%Author: Andrew Beshay, zxxxxxx
%Function for the estimation distance between OOIs.
function r = ExtractOOIs(X, Y, intensities, Selection)
    %Distances
    r.Distance = sqrt(diff(X).^2 + diff(Y).^2);
    %Indexes of each cluster
    r.Clusters = [0 find(r.Distance > Selection.Threshold)'];
    array = length(r.Clusters)-1;
    r.N = 0;
    r.Identity = [];
    r.Id = [0];
    r.Colours   = zeros(1,array);
    r.Centers   = zeros(2,array);
    r.Diameters = zeros(1,array);
    
    for j=1:array
        switch Selection.Case
            case 0
                %DIST method
                XY = [X(r.Clusters(j)+1:r.Clusters(j+1)) Y(r.Clusters(j)+1:r.Clusters(j+1))];
                temp = [XY(1,1) XY(1,2); XY(end,1) XY(end, 2)];
                r.Diameters(j) = sqrt(range(XY(:,1))^2+range(XY(:,2))^2);%pdist(temp);
                r.Centers(:,j) = [mymean(XY(:,1)), mymean(XY(:,2))];%[mean(XY(:, 1)), mean(XY(:,2))];
                %Manually made my own mean function to improve processing
                %time, as well as resorting to a manual dist calculation
                %instead of using pdist.
            case 2
                 %CIRCLEFIT method (Author David Pratt)
                XY = [X(r.Clusters(j)+1:r.Clusters(j+1)) Y(r.Clusters(j)+1:r.Clusters(j+1))];
                if length(XY) >= 3
                    par = CircleFitByPratt(XY);
                    r.Centers(:,j) = [par(1), par(2)];
                    r.Diameters(j) = 2*par(3);
                else 
                    r.Diameters(j) = -1;
                end
            case 1 
                %CIRCLEFIT method (Author Tutorail)
                XY = [X(r.Clusters(j)+1:r.Clusters(j+1)) Y(r.Clusters(j)+1:r.Clusters(j+1))];
                if length(XY) >= 3
                    [temp(1), temp(2), temp(3)] = circlefit(XY(:,1),XY(:,2));   
                    r.Centers(:,j) = [temp(1), temp(2)];
                    r.Diameters(j) = 2*temp(3);
                else 
                    r.Diameters(j) = -1;
                end
            case 3
                %Range method
%                 XY = [X(r.Clusters(j)+1:r.Clusters(j+1)) Y(r.Clusters(j)+1:r.Clusters(j+1))];
%                 temp_x = XY(:,1)';
%                 temp_y = XY(:,2)';
%                 r.Diameters(j) = sqrt(abs(range(temp_x)^2+range(temp_y)^2));
%                 r.Centers(:,j) = [mean(XY(:, 1)), mean(XY(:,2))];
                object_X = X(r.Clusters(j)+1:r.Clusters(j+1));
                object_Y = Y(r.Clusters(j)+1:r.Clusters(j+1));
                r.Centers(:,j) = [mean(object_X),mean(object_Y)];
                r.Diameters(j) = norm(range(object_X),range(object_Y));
                
        end
        if (r.Diameters(j) > 0.2 || r.Diameters(j) < 0.05)
            r.Diameters(j) = -1;
        end    
        if max(intensities(r.Clusters(j)+1:r.Clusters(j+1)))
            r.Colours(j) = 1;
            r.Id = r.Id + 1;
            r.Id = [0 r.Id];
        else
            r.Colours(j) = 0;
        end
    end
    
    r.Centers(:,r.Diameters < 0) = []; 
    r.Colours(r.Diameters < 0) = [];
    r.Diameters(r.Diameters < 0) = [];
    r.N = length(r.Diameters);
%     r.Centers2 = 
    r.Identity = length(r.Centers(1,find(r.Colours)));
    r.Id = nonzeros(r.Id); r.Id = r.Id';
    r.tag = [1:r.N];  

return;
end

function [um] = mymean(x)
    t = length(x);
    sum = 0;
    for i=1:t
        sum = sum + x(i);
    end
    um = sum./t;
    
return;
end

function [xc,yc,rad] = circlefit(X,Y)
    %From the tutors in the tutorial board
    A = -2*[X(1:end-1) - X(2:end), Y(1:end-1)-Y(2:end)];
    B = [X(2:end).^2 - X(1:end-1).^2 + Y(2:end).^2 - Y(1:end-1).^2];
    center = A\B;
    xc = center(1);
    yc = center(2);
    rad = mean(sqrt((X - center(1)).^2 + (Y - center(2)).^2));
return;
end