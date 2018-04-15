## Project: Kidnapped Vehicle Project [![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

Overview
---
This project implements the Kidnappped Vehicle Project for Udacity Self Driving Nano Degree Term 2. It localizes the position and bearing of a vehicle using an initial noisy GPS estimate, as well as a series of sensor observations. It compares the sensor observations to a list of stored landmarks, and uses particle filters to zero in on the position of the vehicle inside a global coordinate system.

Code changes
---
For this project I made the following changes to the starter code.

1. __*particle_filter.cpp*__

1. Initialized a vector of 100 particles drawn from a Gaussian distribution centred around the initial GPS position, and with a standard-deviation corresponding to the sensor tolerances. Each particle is assigned a unique identifier for bookkeeping purposes, and the weights of all particles is set to equal 0.01. I experimented with different numbers of particles (10, 100, 1000, 2000) and found that using 100 particles gives a good tradeoff of execution speed versus localization error.
2. I implemented the prediction method of the ParticleFilter class. This method takes the position and bearing of each of the particles in my particle filter, and updates its position and bearing based on noisy control inputs. The noise is modeled as a Gaussian distribution centred around the particles predicted position (assuming noiseless control), with a standard deviation corresponding to the controls tolerances.
3. I implemented the dataAssociation method which takes in a vector of landmarks (from a stored map) which are within sensor range of a particles position hypothesis, and associates the closest sensor observation to that landmark. This association will be used in assigning weights to each particle within the updateWeights method (discussed next).
4. I implemented the updateWeights method by transforming the observations from vehicle coordinates to map coordinates for each particle, associating the closest prediction to each landmark observation within sensor range of each particles position and estimating a weight for the particle based on a multivariate Gaussian distribution. I also normalized the weights across all particles so the weights can be interpreted as a probability distribution function.
5. For each particle, I also track the landmark IDs of its associations as well as the closest sensor measurment corresponding to each of the associated landmarks. This is used by the simulator to depict how closely the best-match particle tracks the location of the closest landmarks.
6. I implemented the resample method which draws particles from my initial set of particles with a likelihood proportional to the particles weight.

Performance
---
The video below shows the performance of the particle filter. As shown in the video, the performance of the particle filter is better than the thresholds stated in the Project rubric for localization error in position and yaw. Also the simulation completes in less than half the time required by the project.

[![Particle Filter Simulation](https://github.com/calvinhobbes119/Kidnappped-Vehicle-Project/blob/master/Untitled.png)](https://youtu.be/DiNmNsk03V0)

I also implemented a visualization of the particle filter points before and after resampling. The particle filter points before resampling are show in blue, and the after-resampling points are shown in red.

[![Particle Filter Points](https://github.com/calvinhobbes119/Kidnappped-Vehicle-Project/blob/master/Untitled1.png)](https://youtu.be/8oSff9zLmrU)
