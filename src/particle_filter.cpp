/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <typeinfo>
#include <limits>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
  
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  num_particles = 100;  // TODO: Set the number of particles
  std::default_random_engine gen; //Random generator
  
  // Create a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
    
  for(int i=0;i<num_particles;i++){
    //Sample from the normal distributions and assign to ParticleFilter instance variables
    Particle p;
    
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    
    //Initialize all weights to 1
    p.weight = 1.0;
    
    //Add the particle filter to the vector of particles
    particles.push_back(p);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen; //Random generator
  double const1,const2;

  normal_distribution<double> dist_x(0.0 , std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);  
  
  //Using the bicycle motion model calcluate the new positions and heading
  for(int i=0;i<num_particles;i++){
     
    if(fabs(yaw_rate) < 0.0001){  //Avoid division by zero
       particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
       particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
    } 
    else{
      const1 = velocity/yaw_rate;
      const2 = yaw_rate * delta_t;
      particles[i].x = particles[i].x + const1*(sin(particles[i].theta + const2) - sin(particles[i].theta));
      particles[i].y = particles[i].y + const1*(cos(particles[i].theta) - cos(particles[i].theta + const2));
      particles[i].theta = particles[i].theta + const2;
    }
    
    //Add gaussian noise to the predictions
    particles[i].x +=  dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  double distance,min_dist;
  for(int j=0;j<observations.size();j++){
     //Perform data association by finding closest landmark for each observation
     min_dist = std::numeric_limits<double>::max();
     for(int k=0;k<predicted.size();k++){
        distance = dist(observations[j].x, observations[j].y, predicted[k].x, predicted[k].y);
        if( distance < min_dist){
           min_dist = distance;
           observations[j].id = predicted[k].id;
        }   
     }
  }     
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
   //For each partcile find the predicted measurements for landmarks within the sensor range
   vector<LandmarkObs> predictions;
   vector<LandmarkObs> observations_transformed;
   double x_map,y_map,distance;
  
   for(int i=0;i<num_particles;i++){
      predictions.clear();
      for(int j=0;j<map_landmarks.landmark_list.size();j++){
        LandmarkObs obj1;
        obj1.id = map_landmarks.landmark_list[j].id_i;
        obj1.x = map_landmarks.landmark_list[j].x_f;
        obj1.y = map_landmarks.landmark_list[j].y_f;
        distance = dist(particles[i].x, particles[i].y,obj1.x,obj1.y);
        if(distance <= sensor_range)
           predictions.push_back(obj1);   
      }
      //For the predicted measurements and particles, complete data association for each particle
      //First the observations needs to be transformed from vehicle co-ordinates to map co-ordinates
      observations_transformed.clear();
      for(int k=0;k<observations.size();k++){
         x_map = particles[i].x + (cos(particles[i].theta) * observations[k].x) - (sin(particles[i].theta) * observations[k].y);
         y_map = particles[i].y + (sin(particles[i].theta) * observations[k].x) + (cos(particles[i].theta) * observations[k].y);
         observations_transformed.push_back(LandmarkObs{observations[k].id,x_map,y_map});
      }    

     dataAssociation(predictions, observations_transformed);

     //Calculate the weights
     double mu_x, mu_y;
     //Re-initialize weight
     particles[i].weight = 1.0;
     for(int j=0;j<observations_transformed.size();j++){
        for(int k=0;k<predictions.size();k++){
           if(predictions[k].id == observations_transformed[j].id){
              mu_x = predictions[k].x;
              mu_y = predictions[k].y;
              break;
           }
        }
      // Calculate weight using multivariate Gaussian 
      particles[i].weight *= multiv_prob(std_landmark[0], std_landmark[1], observations_transformed[j].x, observations_transformed[j].y, mu_x, mu_y); 
     }     
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   std::default_random_engine gen; //Random generator
   vector<Particle> resampled_particles;
   vector<double> particle_weights;
   resampled_particles.clear();
   particle_weights.clear();
   
  for(int i=0;i<num_particles;i++){
     particle_weights.push_back(particles[i].weight);
  }  
   
  /*std::discrete_distribution produces random integers on the interval [0, n), where the probability 
  of each individual integer i is defined as w i/S, that is the weight of the ith integer divided by the sum of all n weights.*/
  // Reference: https://knowledge.udacity.com/questions/256628
  
  std::discrete_distribution<size_t> dist_index(particle_weights.begin(), particle_weights.end());
  for(int i=0;i<num_particles;i++){
     resampled_particles.push_back(particles[dist_index(gen)]);
  }  
  particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
 
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}