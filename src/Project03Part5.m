%Author: Andrew Beshay, Zxxxxxx
%Program: Solution for AAS, T1.2019, Project3.Part5

function Project03Part5
    %% Loading all prject files
    F = dir('*.mat');
    for i=1:length(F)
        load([F(i).name]);
    end
    clear F;
    %% Master Variables to be used throughout the program
    global Master;
    Master.Shutdown = 0;
    Master.Pause = 0;
    Master.Case = 0;
    Master.Threshold = 0.08;
    Master.DisplayLaser = 0;
    %% Initialising Plots
    Handle = IniPlots();
    %% Importing Data
    % IMU
    IMU.Time = double((IMU.times - IMU.times(1)))*0.0001;
    YawRate = -IMU.DATAf(6,:);
    % Velocity
    Velocity = Vel.speeds;
    % Laser Data
    Laser.Time = double((dataL.times - dataL.times(1)))*0.0001;
    Laser.Scans = dataL.Scans;
    Laser.Angles = dataL.angles;
    CurrentScan = 1;
    %% Initialisation
    Gyro = [0 0 pi/2];
    Yaw = zeros(IMU.N, 1);
    X = zeros(IMU.N, 1);
    Y = zeros(IMU.N, 1);
    Yaw(1) = Gyro(3);
    
    Xe = [0; 0; pi/2; 0; 0];
    P = zeros(5,5);
    P(4,4) = (2*pi/180)^2;
    
    Xe_History = zeros(5,IMU.N);
    
    stdDevGyro = 1.4*pi/180;
    stdDevSpeed = 0.4;
    P(5,5) = stdDevSpeed^2;
    sdev_rangeMeasurement = 0.16;
    sdev_bearingMeasurement = 1.1*pi/180;    
    Q = diag([(0.005)^2,(0.005)^2,(0.01*pi/180)^2,(1*pi/180*1/600)^2,(0.005*1.2)^2]);
    Pu = diag(stdDevGyro^2);
    %% Mean bias (Project 2 method)
    bias = mean(YawRate(1:find(IMU.Time == 17)));
    %% Global Coordinate Frame Initialisation
    % Global coordinate frame is the one alligned with the platform at time
    % 0
    [R, I] = ExtractLaserData(Laser.Scans(:,1));
    Laser.X = -cos(Laser.Angles).*R;
    Laser.Y =  sin(Laser.Angles).*R + 0.46;
    
    GOOIs = ExtractOOIs(Laser.X, Laser.Y, I, Master);
    GOOIs.Centers = GOOIs.Centers(:,find(GOOIs.Colours));
    GOOIs.Centers = Transform(Yaw(1)-pi/2, X(1), Y(1),...
        GOOIs.Centers(1,:), GOOIs.Centers(2,:));
    set(Handle.GOOIs, 'xdata', GOOIs.Centers(1,:), 'ydata', GOOIs.Centers(2,:));
    for i=1:GOOIs.Identity
        set(Handle.G_Labels(i), 'Position', [GOOIs.Centers(1,i)-0.5 GOOIs.Centers(2,i)+0.02],...-0.5],...
            'String', ['\bf',num2str(i)]);
    end
    %% Dynamic Local Coordiante Frame Calculation    
    pause(1);
    plotspeed = 50;
    for i=1:IMU.N-1
        while Master.Pause, pause(0.2); continue; end
        if (Master.Shutdown), break; end

        Dt = IMU.Time(i+1) - IMU.Time(i);
        Yaw(i+1) = Yaw(i) + Dt*(YawRate(i)-bias);
        X(i+1) = X(i) + Dt*Velocity(i)*cos(Yaw(i));
        Y(i+1) = Y(i) + Dt*Velocity(i)*sin(Yaw(i));
        
        Xe = Xe + Dt*[Xe(5)*cos(Xe(3)); Xe(5)*sin(Xe(3)); YawRate(i)-Xe(4); 0; 0];
        
        J = [1,0,-Dt*Xe(5)*sin(Xe(3)),   0, Dt*cos(Xe(3)); 
             0,1, Dt*Xe(5)*cos(Xe(3)),   0, Dt*sin(Xe(3)); 
             0,0,             1      , -Dt,       0      ; 
             0,0,             0      ,   1,       0      ;
             0,0,             0      ,   0,       1     ];
        Fu = [0; 0; Dt; 0; 0];
        Qu = Fu*Pu*Fu';
        P = J*P*J' + Q + Qu;
        
        if (IMU.Time(i) > Laser.Time(CurrentScan) && IMU.Time(i) <= Laser.Time(end))
            CurrentScan = CurrentScan + 1;
            Alpha = Yaw(i) - pi/2;
            alpha = Xe(3) - pi/2;
            [R, I] = ExtractLaserData(Laser.Scans(:,CurrentScan));
            
            s = sprintf('Laserscan: %d IMUscan: %d Time: %.2fs',CurrentScan, i, IMU.Time(i));
            set(Handle.title, 'String', s);
            
            Laser.X = -cos(Laser.Angles).*R;
            Laser.Y =  sin(Laser.Angles).*R + 0.46;
            % Extract Local Coordiante Data
            LOOIs = ExtractOOIs(Laser.X, Laser.Y, I, Master);
            EKF.Centers = LOOIs.Centers(:,find(LOOIs.Colours));
            LOOIs.Centers = LOOIs.Centers(:, find(LOOIs.Colours));
            LOOIs.Centers = Transform(Alpha, X(i), Y(i), LOOIs.Centers(1,:),...
                LOOIs.Centers(2,:));
            EKFS.Centers = Transform(alpha, Xe(1), Xe(2), EKF.Centers(1,:), EKF.Centers(2,:));
            set(Handle.LOOIs, 'xdata', LOOIs.Centers(1,:), 'ydata', LOOIs.Centers(2,:));
            set(Handle.EOOIs, 'xdata', EKFS.Centers(1,:), 'ydata', EKFS.Centers(2,:));

            MeasuredRanges = []; MeasuredBearings = []; 
            if ~isempty(EKFS.Centers(1,:))
                [indexes, nDetectedLandmarks, str] = DataAssociationEKF(EKFS, GOOIs, Handle.EKF_Labels, Handle.LinesEKF);
                EKFIndex = indexes.EKF;
                GIndex = indexes.GLOB;             
                for k=1:nDetectedLandmarks
                    b = EKFIndex(k);
                    dx = EKFS.Centers(1,b) - Xe(1)-0.46*cos(Xe(3));
                    dy = EKFS.Centers(2,b) - Xe(2)-0.46*sin(Xe(3));
                    MeasuredRanges = [MeasuredRanges, sqrt((dx.*dx + dy.*dy))]; 
                    MeasuredBearings = [MeasuredBearings, atan2(dy,dx) + pi/2 - Xe(3)];
                end
                
                for j=1:nDetectedLandmarks
                    a = GIndex(j);
                    b = EKFIndex(j);

                    eDX = GOOIs.Centers(1,a) - Xe(1)-0.46*cos(Xe(3));
                    eDY = GOOIs.Centers(2,a) - Xe(2)-0.46*sin(Xe(3));
                    eDD = sqrt(eDX*eDX + eDY*eDY);
                    eDB = atan2(eDY,eDX) + pi/2 - Xe(3);

                    H = [-eDX/eDD, -eDY/eDD, 0, 0, 0;
                    eDY/eDD^2, -eDX/eDD^2, -1, 0, 0];
                    ExpectedRange = eDD ;  
                    ExpectedBearing = eDB;

                    z = [MeasuredRanges(j) - ExpectedRange;
                    wrapToPi(MeasuredBearings(j) - ExpectedBearing)];
                    R = 4*diag([(sdev_rangeMeasurement)^2, (sdev_bearingMeasurement)^2]);
                    S = R + H*P*H';
                    iS = inv(S);
                    K = P*H'*iS;
                    Xe= Xe+K*z;
                    P = P-K*H*P;
                end
            
            end

            % Current Laserscan Data
            if Master.DisplayLaser
                LaserData = Transform(Alpha, X(i), Y(i), Laser.X', Laser.Y');
                set(Handle.LaserScan, 'xdata', LaserData(1,:), 'ydata', LaserData(2,:));
            else
                set(Handle.LaserScan, 'xdata', 0, 'ydata', 0);
            end
            
        end

        Xe_History(:,i) = Xe;
        addpoints(Handle.veh,X(i),Y(i));
        addpoints(Handle.ekfveh,Xe_History(1,i),Xe_History(2,i));
        addpoints(Handle.bias,IMU.Time(i),rad2deg(Xe_History(4,i)));
        addpoints(Handle.velekf,IMU.Time(i),Xe_History(5,i));
        addpoints(Handle.vel,IMU.Time(i),double(Velocity(i)));

        if mod(i,plotspeed) == 0
            drawnow;
        end
    end
    drawnow
end
function [indexes, landmarks, str] = DataAssociationEKF(L, G, Lable, Line)
    indexes.EKF = zeros(1, length(L.Centers(1,:)));
    indexes.GLOB = zeros(1, length(L.Centers(1,:)));
    landmarks = 0; Lines = zeros(4, 5); str = "Error EKF"; 
    for j=1:G.Identity
        distx = L.Centers(1,:) - G.Centers(1,j);
        disty = L.Centers(2,:) - G.Centers(2,j);
        dist = sqrt(distx.^2 + disty.^2);
        [mindist, id] = min(dist);
        if mindist < 0.4
%             set(Lable(j),'Position',[G.Centers(1,j) + 0.1,...
%                 G.Centers(2,j) - 0.25],'String',...
%                 ['\bf', '\fontsize{15}\epsilon\fontsize{10}: ',...
%                 num2str(mindist,'%.4f m #'),num2str(j)]);
                set(Lable(j),'Position',[G.Centers(1,j) + 0.1,...
                G.Centers(2,j) - 0.25],'String',...
                ['\fontsize{10}','#',num2str(j)]);
            str = str + newline + ['#',num2str(j),' ',num2str(mindist,'%.4f m')];
            landmarks = landmarks + 1;
            indexes.GLOB(landmarks) = j;
            indexes.EKF(landmarks) = id;     
            Lines(1:2,j) = [G.Centers(1,j); L.Centers(1,id)];%G.Centers(1,j) + distx(id)];
            Lines(3:4,j) = [G.Centers(2,j); L.Centers(2,id)];%G.Centers(2,j) + distx(id)];
            set(Line(j),'xdata',[Lines(1,j)' Lines(2,j)],...
                'ydata',[Lines(3,j)' Lines(4,j)']);
        else
            set(Lable(j),'Position',[0 0],'String','');
            set(Line(j),'xdata',0,'ydata',0);
        end
    end
end
%% ----------------------------------------------------
% Formula/Equations given by Demonstrators in Tutorials
function result = Transform(alpha, Xr, Yr, X, Y)
    result = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]*...
        [X; Y];
    result = result + [Xr; Yr];
end
%% ----------------------------------------------------
function [Ranges, Intensities] = ExtractLaserData(scans)
    MaskLow13Bits   = uint16(2^13-1);
    MaskHigh13Bits  = bitshift(uint16(7), 13);
    RangesA         = bitand(scans, MaskLow13Bits);
    Ranges          = 0.01*double(RangesA);
    Intensities     = bitand(scans, MaskHigh13Bits);  
end
%% ----------------------------------------------------
function m = IniPlots()
    global Master;
    
    m.f1 = figure(1); clf(); hold on;
    set(gca,'color',[0.23, 0.23, 0.23]);
    ax = gca;
    ax.GridColor = [0.93, 0.93, 0.93];
    m.title = title('');
    m.ekfveh = animatedline('color','g');
    m.veh = animatedline('color','c');
    m.LaserScan = plot(0,0,'w.'); 
    m.GOOIs = plot(0,0,'o','color',[1.0 0.5 0.6],'LineWidth',1.1,'MarkerSize',6);
    m.LOOIs = plot(0,0,'x','color',[0.2 1.0 0.0],'LineWidth',1.0,'MarkerSize',6); 
    m.EOOIs = plot(0,0,'+','color',[1.0 1.0 0.0],'LineWidth',1.0,'MarkerSize',6); 
    
    uicontrol('Style','togglebutton','String','Pause','Position',[10,2,80,20],'Callback',{@MyCallBackA,1});
    uicontrol('Style','pushbutton','String','End Now','Position',[90,2,80,20],'Callback',{@MyCallBackA,2});
    uicontrol('Style','radiobutton','String','Show all','Position',[360,2,80,20],'Callback',{@MyCallBackA,3});
    uicontrol('Style','popupmenu','String',{'Distance','CircleFitByTute', 'CircleFitByPratt'},...
        'Position',[470,2,80,20],'Callback',{@MyCallBackA,4});
    set( findall( m.f1, '-property', 'Units' ), 'Units', 'Normalized' )
    
    m.G_Labels = text(zeros(1,5),zeros(1,5),'','Color',[0.93 0.93 0.93]); 
    m.L_Labels = text(zeros(1,5),zeros(1,5),'','Color',[0.5, 0.85, 1.0]); 
    m.EKF_Labels = text(zeros(1,5),zeros(1,5),'','Color',[0.5, 0.85, 1.0]);
    m.LinesEKF = gobjects(5,1);
    m.LinesL = gobjects(5,1);
    for i=1:5
       m.LinesEKF(i) = line(0,0,'color',[0.503, 0.4, 1.0]);
       m.LinesL(i) = line(0,0,'color',[0.503, 0.4, 1.0]);
    end

    xlabel('X (meters)'); ylabel('Y (meters)'); grid on;
    m.legend = legend('\color[rgb]{0.93, 0.93, 0.93}Veh Dr',...
        '\color[rgb]{0.93, 0.93, 0.93}Veh EKF',...
        '\color[rgb]{0.93, 0.93, 0.93}LaserScan',...
        '\color[rgb]{0.93, 0.93, 0.93}Landmarks',...
        '\color[rgb]{0.93, 0.93, 0.93}OOIs',...
        '\color[rgb]{0.93, 0.93, 0.93}EKF',...
        '\color[rgb]{0.93, 0.93, 0.93}error','AutoUpdate','off');

    m.legend.Color = [0.38, 0.38, 0.38];
    m.legend.Title.Color = [0.93, 0.93, 0.93];
    axis([-5, 4, 0, 7]);
    zoom on;
     
    m.f2 = figure(2); clf(); hold on; grid on;
    m.sb1 = subplot(2,1,1); grid on;
    m.title2 = 'Bias error';
    xlabel('Time (seconds)'); ylabel('Bias (degrees)');
    set(gca,'color',[0.23, 0.23, 0.23]);
    ax = gca;
    ax.GridColor = [0.93, 0.93, 0.93];
    m.bias = animatedline('color',[0.5, 0.85, 1.0],'MaximumNumPoints',2500);
    axis tight;
    
%     m.f3 = figure(3); clf(); hold on; grid on;
    m.sb2 = subplot(2,1,2); grid on;
    m.title3 = 'Velocity Estimation';
    xlabel('Time (seconds)'); ylabel('Velocity (m/s)');
    set(gca,'color',[0.23, 0.23, 0.23]);
    ax = gca;
    ax.GridColor = [0.93, 0.93, 0.93];
    m.vel = animatedline('color',[0.8, 0.85, 0.4]);%,'MaximumNumPoints',2500);
    m.velekf = animatedline('color',[0.5, 0.85, 1.0]);%,'MaximumNumPoints',2500);
    axis tight;
end
%% ----------------------------------------------------
function MyCallBackA(src,~,x)   
    global Master;
    f1 = figure(1);

    switch x
        case 1
            Master.Pause = ~Master.Pause; 
            if Master.Pause
                src.String = 'Continue';
            else
                src.String = 'Pause';
            end
        case 2
            Master.Shutdown = ~Master.Shutdown;
            disp('Ending the program');
        case 3
            Master.DisplayLaser = ~Master.DisplayLaser;
            if Master.DisplayLaser
                f1.CurrentAxes.XLim = [-10,10];
                f1.CurrentAxes.YLim = [-10,10];
                disp("Showing complete laser scan");
            else
                f1.CurrentAxes.XLim = [-5,4];
                f1.CurrentAxes.YLim = [0,7];
                disp("Showing only OOI's");
            end
        case 4
            val = src.Value;
            switch val
                case 1
                    Master.ToggleMethod = 0;
                case 2
                    Master.ToggleMethod = 1;
                case 3
                    Master.ToggleMethod = 2;
                case 4
                    Master.ToggleMethod = 3;
            end
        case 5
            Master.Threshold = src.Value;
            set(Master.Textbox, 'String', num2str(Master.Threshold));
    end
end