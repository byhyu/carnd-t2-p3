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

#include "particle_filter.h"

using namespace std;

#define DEBUG

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    #ifdef DEBUG
    cout << "======== Init Particles ====== " << endl;
    #endif

    num_particles = 100;
    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

// Init weigths with 1.0
    weights.resize(num_particles);
    fill(weights.begin(), weights.end(), 1.0);

    // init particles
    for (int i = 0; i <num_particles; i++) {
        Particle p;
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = weights[i];
        particles.push_back(p);
    }
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

#ifdef DEBUG
    cout << "======== Prediction ====== " << endl;
#endif

    default_random_engine gen;
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
    normal_distribution<double> dist_theta(0.0, std_pos[2]);

    for (int i = 0; i < num_particles; i++) {

#ifdef DEBUG
        printf("Particle %4d: (%f, %f, %f)->", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
#endif

        double travel_dist = 0; // either yaw dist or forward dist
        if (fabs(yaw_rate) <= 0.00001) {
            travel_dist = velocity * delta_t;  // distance
            particles[i].x += travel_dist * cos(particles[i].theta) + dist_x(gen);
            particles[i].y += travel_dist * sin(particles[i].theta) + dist_y(gen);
            particles[i].theta += dist_theta(gen);
        } else {
            double v_scale = velocity / yaw_rate;
            travel_dist = yaw_rate * delta_t;  // rotational distance
            particles[i].x += v_scale * (sin(particles[i].theta + travel_dist) - sin(particles[i].theta)) + dist_x(gen);
            particles[i].y += v_scale * (cos(particles[i].theta) - cos(particles[i].theta + travel_dist)) + dist_y(gen);
            particles[i].theta += travel_dist + dist_theta(gen);
        }
        // particles[i].theta = fmod(particles[i].theta, 2 * M_PI);  //TODO: Normalize or not?

        #ifdef DEBUG
        printf("(%f, %f, %f)\n", particles[i].x, particles[i].y, particles[i].theta);
        #endif

    }
#ifdef DEBUG
    cout << "finished prediction" << endl;
#endif

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

#ifdef DEBUG
    cout << "======= data association =======";
#endif

    for (int i = 0; i < observations.size(); i++) {
        double curr_dist = 0;
        double closest_dist = numeric_limits<double>::max();
        for (int j = 0; j < predicted.size(); j++) {
            curr_dist = fabs(dist(predicted[j].x,
                                  predicted[j].y,
                                  observations[i].x,
                                  observations[i].y));
            if (curr_dist < closest_dist) {
                observations[i].id = predicted[j].id;
                closest_dist = curr_dist;
            }
        }

    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html

#ifdef DEBUG
    cout << "======= Update Weights =======" << endl;

#endif

    for (int i = 0; i < num_particles; i++) {
        vector<LandmarkObs> predicted_locs;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            if (abs(dist(map_landmarks.landmark_list[j].x_f,
                         map_landmarks.landmark_list[j].y_f,
                         particles[i].x,
                         particles[i].y)) <= sensor_range) {
                LandmarkObs l;
                l.id = map_landmarks.landmark_list[j].id_i;
                l.x = map_landmarks.landmark_list[j].x_f;
                l.y = map_landmarks.landmark_list[j].y_f;
                predicted_locs.push_back(l);

            }
        }

        // transform sensor measurements to MAP coord system
        vector<LandmarkObs> transformed_obs;
        for (int j = 0; j < observations.size(); ++j) {
            LandmarkObs l;
            l.id = observations[j].id;
            l.x = particles[i].x + cos(particles[i].theta) * observations[j].x -
                  sin(particles[i].theta) * observations[j].y;
            l.y = particles[i].y + sin(particles[i].theta) * observations[j].x +
                  cos(particles[i].theta) * observations[j].y;
            transformed_obs.push_back(l);
        }

        // find associations between obs and map landmarks
        dataAssociation(predicted_locs, transformed_obs);
        particles[i].associations.clear();
        particles[i].sense_x.clear();
        particles[i].sense_y.clear();
        for (int j = 0; j > transformed_obs.size(); ++j) {
            particles[i].associations.push_back(transformed_obs[j].id);
            particles[i].sense_x.push_back(transformed_obs[j].x);
            particles[i].sense_y.push_back(transformed_obs[j].y);
        }

        // update weights based on associations
        // calculate normalization term
        double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
        particles[i].weight = 1;

        for (int j = 0; j < transformed_obs.size(); ++j) {

#ifdef  DEBUG
            cout << "update weights: ++j" << j << endl;
#endif

            double x_obs = transformed_obs[j].x;
            double y_obs = transformed_obs[j].y;
            double mu_x = map_landmarks.landmark_list[transformed_obs[j].id - 1].x_f;
            double mu_y = map_landmarks.landmark_list[transformed_obs[j].id - 1].y_f;
            // calculate exponent
            double exponent = (pow(x_obs - mu_x, 2) / (2 * pow(std_landmark[0], 2))) +
                              (pow(y_obs - mu_y, 2) / (2 * pow(std_landmark[1], 2)));
            // calculate weight using normalization terms and exponent
            particles[i].weight *= gauss_norm * exp(-exponent);

        }

#ifdef DEBUG
        cout << "Weights updated" << endl;
#endif

    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

#ifdef DEBUG
    cout << "====== resample =======" << endl;
#endif

    weights.clear();
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }

    // create sampler
    default_random_engine gen;
    discrete_distribution<> d(weights.begin(), weights.end());

    // resample
    vector<Particle> resampled_particles;
    for (int i = 0; i < num_particles; i++) {
        Particle p = particles[d(gen)];
        p.id = i;
        resampled_particles.push_back(p);

#ifdef DEBUG_OUTPUT
        printf("Particle %4d: (%f, %f, %f, %f)->(%f, %f, %f)\n", particles[i].id, particles[i].x, particles[i].y, particles[i].theta, weights[i], p.x, p.y, p.theta);
#endif
    }
    particles.swap(resampled_particles);
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
