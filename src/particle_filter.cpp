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
    for (int i = 0; i <num_particles; ++i) {
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
    normal_distribution<double> dist_x(x, std_pos[0]);
    normal_distribution<double> dist_y(y, std_pos[1]);
    normal_distribution<double> dist_theta(theta, std_pos[2]);

    for (int i = 0; i < num_particles; ++i) {

#ifdef DEBUG
        printf("Particle %4d: (%f, %f, %f)->", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
#endif

        double travel_dist = 0; // either yaw dist or forward dist
        if (abs(yaw_rate) <= 1e-10) {
            travel_dist = velocity * delta_t;  // distance
            particles[i].x += travel_dist * cos(particles[i].theta) + dist_x(gen);
            particles[i].y += travel_dist * sin(particles[i].theta) + dist_y(gen);
            particles[i].theta += dist_theta(gen);
        } else {
            double v_scale = velocity / yaw_rate;
            travel_dist = yaw_rate * delta_t;  // rotational distance
            particles[i].x += v_scale * (sin(particles[i].theta + travel_dist) - sin(particles[i].theta)) + nd_x(gen);
            particles[i].y += v_scale * (cos(particles[i].theta) - cos(particles[i].theta + travel_dist)) + nd_y(gen);
            particles[i].theta += travel_dist + dist_theta(gen);
        }
        particles[i].theta = fmod(particles[i].theta, 2 * M_PI);  //TODO: Normalize or not?

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
            curr_dist = abs(dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y));
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
    cout << "======= Update Weights ======="<<endl;
#endif



}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
