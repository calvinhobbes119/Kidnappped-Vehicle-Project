/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1/num_particles. 
	// Add random Gaussian noise to each particle.
    int i; //index for particle
    std::default_random_engine gen;
    std::normal_distribution<double> x_distribution(x,std[0]), y_distribution(y,std[1]), theta_distribution(theta,std[2]);
    Particle particle_template; // template used for particle creation
    double x_noise, y_noise, theta_noise; //

    num_particles = 100; // Experimentally found that 100 particles 
                         // gives good tradeoff of speed vs accuracy
    particle_template.id = 0;
    particle_template.x = 0;
    particle_template.y = 0;
    particle_template.theta = 0;
    

    for (i=0; i<num_particles;i++)
    {
        Particle new_particle = particle_template;
        
        // Set each particles position and bearing from a
        // Gaussian distribution centred around the GPS measurement
        new_particle.x = x_distribution(gen);
        new_particle.y = y_distribution(gen);
        new_particle.theta = theta_distribution(gen);
        new_particle.id = i; // Assign a unique identifier for each particle
        particles.push_back(new_particle);
        weights.push_back(1.0/num_particles);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    int i; // Index for particle
    std::default_random_engine gen;
    double x_noise, y_noise, theta_noise;
    Particle my_particle;

    for (i=0; i<num_particles;i++)
    {
        my_particle = particles[i];
        if (yaw_rate != 0)
        {
            // Predict the particle's position assuming non-zero yaw rate.
            my_particle.x += velocity/yaw_rate * (sin(my_particle.theta + yaw_rate * delta_t) - sin(my_particle.theta));
            my_particle.y += velocity/yaw_rate * (-cos(my_particle.theta + yaw_rate * delta_t) + cos(my_particle.theta));
            my_particle.theta += yaw_rate * delta_t; 
        }
        else
        {
            // Predict the particle's position assuming zero yaw rate. So yaw parameter is not updated here.
            my_particle.x += velocity * (cos(my_particle.theta)) * delta_t;
            my_particle.y += velocity * (sin(my_particle.theta)) * delta_t;
        }
        std::normal_distribution<double> x_distribution(my_particle.x,std_pos[0]), y_distribution(my_particle.y,std_pos[1]), theta_distribution(my_particle.theta,std_pos[2]);
        // Update the particle's position by adding Gaussian noise for motion uncertainity.
        my_particle.x = x_distribution(gen);
        my_particle.y = y_distribution(gen);
        my_particle.theta = theta_distribution(gen);
        particles[i] = my_particle;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	// observed measurement to this particular landmark.
    int i,j; //Index for predicted landmark position and actual landmark observations
    LandmarkObs observation, prediction;
    double min_error;
    double error;
    
    for (i=0; i<observations.size();i++)
    {
        min_error = -1.0;
        observation = observations[i];

        for (j=0; j<predicted.size();j++)
        {
            prediction = predicted[j];
            error = (observation.x - prediction.x)*(observation.x - prediction.x) + 
                    (observation.y - prediction.y)*(observation.y - prediction.y);
            if ((min_error == -1) ||
                (min_error > error))
            {
                // Assign the predicted measurement which is closest match
                // to this observation and update the error.
                observations[i].id = j;
                min_error = error;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution.
	// NOTE: The observations are given in the VEHICLE'S coordinate system. The particles are located
	//   according to the MAP'S coordinate system. We will need to transform between the two systems.
    int i, j; // Indices for particles and observations.
    double x_obs, y_obs, mu_x, mu_y, sum = 0;
    cout << endl;
    for (i=0; i<num_particles; i++)
    {
        std::vector<LandmarkObs> observations_map_coordinates, map_landmarks_within_sensor_range;
        particles[i].sense_x.clear();
        particles[i].sense_y.clear();
        particles[i].associations.clear();

        for (j=0; j < observations.size(); j++)
        {
            LandmarkObs observation_map_coordinates;
            // Translate each observation to the map's coordinate system
            observation_map_coordinates.id = j;
            observation_map_coordinates.x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
            observation_map_coordinates.y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
            observations_map_coordinates.push_back(observation_map_coordinates);
        }
        for (j=0; j < map_landmarks.landmark_list.size(); j++)
        {
            if (dist(particles[i].x, particles[i].y, (double)map_landmarks.landmark_list[j].x_f, (double)map_landmarks.landmark_list[j].y_f) <= sensor_range)
            {
                // This map landmark is within sensor range for the assumed position of the vehicle
                // as represented by this particle. Add it to the list of map landmarks which should be
                // detectable by the vehicles sensors. We will correlate this list with the list of
                // actual observations next.
                LandmarkObs map_landmark;
                map_landmark.id = map_landmarks.landmark_list[j].id_i;
                map_landmark.x = (double)map_landmarks.landmark_list[j].x_f;
                map_landmark.y = (double)map_landmarks.landmark_list[j].y_f;
                map_landmarks_within_sensor_range.push_back(map_landmark);
                // Store the ID associated with this landmark
                particles[i].associations.push_back(map_landmark.id);
            }
        }
        // Find the closest predicted measurement to each landmark within sensor range
        dataAssociation(observations_map_coordinates, map_landmarks_within_sensor_range);
        // Initialize the weights for each particle to 1.0
        weights[i] = 1.0;

        for (j=0; j < map_landmarks_within_sensor_range.size(); j++)
        {
            x_obs = (double)observations_map_coordinates[map_landmarks_within_sensor_range[j].id].x;
            y_obs = (double)observations_map_coordinates[map_landmarks_within_sensor_range[j].id].y;
            // Store the X/Y sensor observation associated with this landmark
            particles[i].sense_x.push_back((double)x_obs);
            particles[i].sense_y.push_back((double)y_obs);
            mu_x  = map_landmarks_within_sensor_range[j].x;
            mu_y  = map_landmarks_within_sensor_range[j].y;
            double temp;
            temp = (1/(2 * M_PI * std_landmark[0] * std_landmark[1])) * exp(-1 *(((x_obs - mu_x)*(x_obs - mu_x))/(2 * std_landmark[0]*std_landmark[0]) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * std_landmark[1]*std_landmark[1])));
            // Probability that the observed measurement was actually from this landmark assuming this particle
            // represents the true position of the vehicle, and sensor noise is Gaussian distributed around the
            // actual position of the landmark.
            weights[i] *= temp;
        }
        sum += weights[i];
    }
    // Normalize the weights so they represent a probability distribution
    for (i=0; i<num_particles; i++)
    {
        weights[i] /= sum;
    }
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
    std::discrete_distribution<> d(weights.begin(), weights.end());
    // Setup the random bits
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<Particle> resampled_particles;
    std::vector<double> resampled_weights;
    cout << endl;
    for(int n=0; n<weights.size(); ++n)
    {
        // Generate a random index between 0...num_particles-1
        // based on the probability distribution represented by
        // the weights vector.
        int t = d(gen);
        resampled_particles.push_back(particles[t]);
        resampled_weights.push_back(weights[t]);
    }

    // Update the particles list, weights and IDs based on
    // the resampled list.
    for(int n=0; n<weights.size(); ++n)
    {
        particles[n] = resampled_particles[n];
        particles[n].id = n;
        particles[n].weight = resampled_weights[n];
    }
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
