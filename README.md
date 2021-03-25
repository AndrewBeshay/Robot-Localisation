# Robot-Localisation

## Part 1) 
Modify the example “DemoEKF_2019.m” for   
1. Processing bearing (angle) observations (in addition to the range ones which are already implemented in the 
example.)  
2. Using a proper Q matrix, based on the assumed noise in the inputs of the process model.  


Note: for solving item (1) you need to also modify certain parts of the simulator components.  
## Part 2)   
Simulate the existence of bias in the gyroscope’s measurements. Extend the EKF solution in (1) for estimating the 
gyroscope’s bias, as explained in lecture notes “[AAS2019]EKF_Localizer_estimating_bias.pdf”.   


Test your implementation by simulating cases having a bias= -1 degree/second, and having a bias =+1 degree/second.  
## Part 3) 
Use (1) for implementing the EKF based on the data and code used in Project2.Part4.  
This program is based on the data used in Project 2. It is OFF-LINE but using real data, in the same way we solved 
Project2. You are requested to adapt your solution for Project2.Part4, adding the EKF component. The estimates should 
be more accurate than those obtained in P2.4, which were based on simple dead-reckoning.  


Assume the following realistic conditions:  
* Noise in angular rate measurements: standard deviation = 1.4 degrees/second.   
* Noise in speed sensor: standard deviation = 0.4m/s.  
* Noise in range measurements: standard deviation = 0.16m.  
* Noise in bearing measurements: standard deviation = 1.1 degree  

You need also to consider that the LIDAR sensor is located at the front of the platform. There is a document, in the 
lecture notes, which describes how to treat that matter.  
Note: you must remove the bias which is present in the angular rate measurements (this is not required in part 4)  
## Part 4) 
Include the capability implemented in (2) for extending the solution developed in (3).  
Consider that the gyroscope’s bias can be, in the worst case, limited in the range [-2,+2] degree/second.   


Note: you must not remove the bias which is present in the raw angular rate measurements; because your estimation   
process estimates and remove it on run-time.AAS - Project 3 – S1.2019 – Version 1  
## Part 5) 
Modify (4) for estimating longitudinal velocity. This means those measurements are not provided, and, 
consequently, you need to estimate the longitudinal component of the velocity, in addition to the platform’s 2D pose
and the gyroscope’s bias.  


Recommendation: you may solve it first modifying your solutions in pure simulation, for verifying performance; before 
trying it with the real rata. Assume that the maximum acceleration which the platform can achieve is 1.5m/s2
(in forward 
or in reverse gear)  
