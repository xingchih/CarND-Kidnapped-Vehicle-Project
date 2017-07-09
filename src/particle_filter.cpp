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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    default_random_engine gen;
    num_particles = 100;
    // Create a normal (Gaussian) distributions for x, y and theta
    normal_distribution<double> dist_x    (x,     std[0]);
    normal_distribution<double> dist_y    (y,     std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    Particle particle_tmp;

    for(size_t i = 0; i<num_particles; i++){

        particle_tmp.x      = dist_x    (gen);
        particle_tmp.y      = dist_y    (gen);
        particle_tmp.theta  = dist_theta(gen);
        particle_tmp.weight = 1.0;
        particles.push_back(particle_tmp);
        weights.push_back(1.0);
    }
    is_initialized = true;
    cout<<"initialized"<<endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

    for(size_t i = 0; i<particles.size(); i++){

        double x0     = particles[i].x;
        double y0     = particles[i].y;
        double theta0 = particles[i].theta;

        double thetaf = theta0 + yaw_rate*delta_t;
        double xf     = x0 + velocity/yaw_rate*( sin(thetaf) - sin(theta0));
        double yf     = y0 + velocity/yaw_rate*( cos(theta0) - cos(thetaf));

        // calculate new state
        if (fabs(yaw_rate) < 1e-5) {  
            xf     = x0 + velocity*delta_t*cos(theta0);
            yf     = y0 + velocity*delta_t*sin(theta0);
            thetaf = theta0;
        } 
        else{
            thetaf = theta0 + yaw_rate*delta_t;
            xf     = x0     + velocity/yaw_rate*( sin(thetaf) - sin(theta0));
            yf     = y0     + velocity/yaw_rate*( cos(theta0) - cos(thetaf));
        }

        normal_distribution<double> dist_x    (xf,     std_pos[0]);
        normal_distribution<double> dist_y    (yf,     std_pos[1]);
        normal_distribution<double> dist_theta(thetaf, std_pos[2]);

        particles[i].x      = dist_x    (gen);
        particles[i].y      = dist_y    (gen);
        particles[i].theta  = dist_theta(gen);
    } 
    //cout<<"ParticleFilter::prediction"<<endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    //cout<<"ParticleFilter::dataAssociation in "<<endl;
    double min_dist = 1e9;
    double dist     = 1e9;
    double x_o, y_o, x_p, y_p, dx, dy;
    
    for(size_t j = 0; j<observations.size(); j++){
        min_dist = 1e9;
        x_o = observations[j].x;
        y_o = observations[j].y;

        for(size_t k = 0; k<predicted.size(); k++){
            x_p  = predicted[k].x;
            y_p  = predicted[k].y;
            dx   = x_o - x_p;
            dy   = y_o - y_p;
            dist = sqrt(dx*dx + dy*dy);

            if(dist<min_dist){
                min_dist = dist;
                observations[j].id = predicted[k].id;
            }
        }
    }
    //cout<<"ParticleFilter::dataAssociation out"<<endl;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
    
    //cout<<"ParticleFilter::updateWeights start"<<endl;

    double x, y, dx, dy, c_theta, s_theta, p_xy, dist;
    vector<LandmarkObs> prd_in_map; //predicted landmarks in map coord
    vector<LandmarkObs> obs_in_map; //observed  landmarks in map coord
    LandmarkObs         landmark_tmp; 
    
    for(size_t j = 0; j<particles.size(); j++){
        prd_in_map.clear();
        obs_in_map.clear();

        // in map coord, select landmarks that are within sensor range
        for(size_t k = 0; k<map_landmarks.landmark_list.size(); k++){
            x    = map_landmarks.landmark_list[k].x_f;
            y    = map_landmarks.landmark_list[k].y_f;

            dx   = x - particles[j].x;
            dy   = y - particles[j].y;

            // distance between landmark and particle
            dist = sqrt(dx*dx + dy*dy);
            
            // select if distance less than sensor range
            if(dist < sensor_range){
                landmark_tmp.id = map_landmarks.landmark_list[k].id_i;
                landmark_tmp.x  = map_landmarks.landmark_list[k].x_f;
                landmark_tmp.y  = map_landmarks.landmark_list[k].y_f;
                
                prd_in_map.push_back(landmark_tmp);
            }
        }

        // transform observation into map coord
        for(size_t k = 0; k<observations.size(); k++){
            // transform into map coodinates
            c_theta = cos( particles[j].theta );
            s_theta = sin( particles[j].theta );

            // simple 2D transform
            x = particles[j].x + observations[k].x * c_theta - observations[k].y * s_theta;// observation transformed into map coordinates
            y = particles[j].y + observations[k].x * s_theta + observations[k].y * c_theta;

            landmark_tmp.id = observations[k].id;
            landmark_tmp.x  = x;
            landmark_tmp.y  = y;
            
            obs_in_map.push_back(landmark_tmp);
        }
         
        // nearest neighbour association of predicted and observed landmarks, both in map coord
        dataAssociation(prd_in_map, obs_in_map);
 
        // init particle weight with 1.0
        particles[j].weight = 1.0;

        // loop through observed landmarks, 
        // find the corresponded true landmark in map
        // calculate multi-variate gaussian
        // and assign weight
        for(size_t k = 0; k<obs_in_map.size(); k++){
            for(size_t i = 0; i<prd_in_map.size(); i++){
                if(obs_in_map[k].id == prd_in_map[i].id){
                    dx   = obs_in_map[k].x - prd_in_map[i].x;
                    dy   = obs_in_map[k].y - prd_in_map[i].y;
                    //calculate multi-variate gaussian
                    p_xy = 1/(2*M_PI*std_landmark[0]*std_landmark[1])
                           * exp(-(dx*dx/(2*std_landmark[0]*std_landmark[0]) + dy*dy/(2*std_landmark[1]*std_landmark[1])));
                    //multiply weight
                    particles[j].weight *= p_xy;
                }
            }
        }
        weights[j] = particles[j].weight;
    }

    //cout<<"ParticleFilter::updateWeights end"<<endl;
    
}

void ParticleFilter::resample() {
    
    //cout<<"ParticleFilter::resample 0"<<endl;
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    //cout<<"ParticleFilter::resample start"<<endl;
    random_device           rd;
    mt19937                 gen(rd());
    discrete_distribution<> d(weights.begin(), weights.end());
    vector<Particle>        particles_vec;
    Particle                particle_tmp;

    for(size_t k=0; k<particles.size(); k++) {
        particle_tmp = particles[d(gen)];
        particles_vec.push_back(particle_tmp);
    }
    particles = particles_vec;
    //cout<<"ParticleFilter::resample end"<<endl;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
