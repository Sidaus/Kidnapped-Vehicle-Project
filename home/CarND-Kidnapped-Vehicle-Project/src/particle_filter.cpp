/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

//I am really sorry for making such mistakes, I try my best not to make such mistakes at these.................

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::numeric_limits;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  //std::default_random_engine gen;//the random engine to generate the values
  //Making sure that the filter isn't initialized before hand
  if(is_initialized==true)
    return;
  
  num_particles = 1005;  // TODO: Set the number of particles

  //Now lets create the normal distributions
  normal_distribution<double> distx(x, std[0]); //Here x is representing the GPS coordinates
  //similarly for y and theta
  normal_distribution<double> disty(y, std[1]);
  normal_distribution<double> distTheta(theta, std[2]);
  
  //NOw to smaple from our parameters
  for (int i = 0; i<num_particles; i++)
  {
    Particle particle;
    particle.id = i;
    particle.x = distx(gen);
    particle.y = disty(gen);
    particle.theta = distTheta(gen);
    particle.weight = 1.0;
    
    particles.push_back(particle);//Creating a vector
  }
  
  //The filter is now completely initialized 
  is_initialized = true;
  
}

//Now for the prediction
void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    // Extracting standard deviations
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_Theta = std_pos[2];

  // Creating normal distributions
  normal_distribution<double> distx(0, std_x);
  normal_distribution<double> disty(0, std_y);
  normal_distribution<double> distTheta(0, std_Theta);
  
  //Lets make the prediction now
  for (int i=0; i<num_particles; i++)
  {
    double Theta = particles[i].theta;
    
    //Now the values of x and y will depend on the fact if the yaw rate is changing or remains same......
    //So we can make two cases
    if (fabs(yaw_rate) < 0.0000001)//This will most probably won't be used but still it is a good practise to make all the cases 
    {
      particles[i].x = particles[i].x + velocity*delta_t*Theta;
      particles[i].y = particles[i].y + velocity*delta_t*Theta;
      //yaw remains the same as the yaw rate is constant
    }
    else
    {
      particles[i].x = particles[i].x + velocity / yaw_rate*(sin(yaw_rate*delta_t + Theta)-sin(Theta));
      particles[i].y = particles[i].y + velocity / yaw_rate*(cos(Theta)-cos(yaw_rate*delta_t + Theta));
      //here yaw will change as well
      particles[i].theta = particles[i].theta + yaw_rate*delta_t;
    }
    
    //Now to add the noise
    particles[i].x += distx(gen);
    particles[i].y += disty(gen);
    particles[i].theta += distTheta(gen);
  }
 //Prediction step is now complete
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
//This step is especially important as it tells the code which patricles to associate with the landmarks 
  for (unsigned i = 0; i<observations.size(); i++)
  {
    double MinDist = numeric_limits<double>::max();
    int MapID = -1;
    
    for (unsigned j = 0; j<predicted.size(); j++)
    {
      double XDist = observations[i].x - predicted[j].x;
      double YDist = observations[i].y - predicted[j].y;
      
      double Dist = XDist*XDist + YDist*YDist; //Distance between the two points
      //Now to make sure the minimum dist is selected correctly so that the correct particles are selected
      if (Dist<MinDist)
      {
        MinDist = Dist;
        MapID = predicted[j].id;
      }
    }
    observations[i].id = MapID;
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

  for (int i=0; i<num_particles; i++)
  {
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    
    double range = sensor_range*sensor_range;
    vector<LandmarkObs> Landmarks_IN_Range;
    
    for(unsigned j=0; j<map_landmarks.landmark_list.size(); j++)
    {
      //get id of the map and also x and y coordinates of the landmarks
      float landmark_X = map_landmarks.landmark_list[j].x_f;
      float landmark_Y = map_landmarks.landmark_list[j].y_f;
      int Map_ID = map_landmarks.landmark_list[j].id_i;
      
      double DistLmX = x - landmark_X;
      double DistLmY = y - landmark_Y;
      
      if ((DistLmX*DistLmX+DistLmY*DistLmY) < range)
      {
        Landmarks_IN_Range.push_back(LandmarkObs{Map_ID, landmark_X, landmark_Y});
      }
    }
  
  //Now we have to tranform the coordintes since the particles are in the Vehicle coordinte system and the landmarks are in 
  //the map coordinate system
  
    vector<LandmarkObs> changed_Observations;
    for (unsigned j=0; j<observations.size(); j++)
    {
      double transformed_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double transformed_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      changed_Observations.push_back(LandmarkObs{ observations[j].id, transformed_x, transformed_y });
    }

  //Associate the observations to the landmarks
    dataAssociation(Landmarks_IN_Range,changed_Observations);
    particles[i].weight = 1.0;
    
    for (unsigned int j=0; j<changed_Observations.size(); j++)
    {
      double obs_x,obs_y;
      double ObsX = changed_Observations[j].x;
      double ObsY = changed_Observations[j].y;
      int Obs_ID = changed_Observations[j].id;
      
      for (unsigned k=0; k<Landmarks_IN_Range.size(); k++)
      {
        if(Landmarks_IN_Range[k].id==Obs_ID)
        {
          obs_x = Landmarks_IN_Range[k].x;
          obs_y = Landmarks_IN_Range[k].y;
        }
      }
      //NOw to calculate the weights using the Multi-Variate-Gaussian (this sounds very cool, saying this increased my IQ by 50 points)
      
      double sx = std_landmark[0];
      double sy = std_landmark[1];
      double dx = ObsX-obs_x;
      double dy = ObsY-obs_y;
      //Now for the observed weights
      double obsw = ( 1/(2*M_PI*sx*sy)) * exp( -( pow(dx,2)/(2*pow(sx, 2)) + (pow(dy,2)/(2*pow(sy, 2))) ) );
      if (obsw == 0)
      {
        particles[i].weight *= 0.0000001;
      }
      else
      {
        particles[i].weight *= obsw;
      }
    }
  }
}

//The final step
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  //I found the mistakes in the update weight code but cannot do it here, can you please be a bit more specific.....
  //I am the not able to find the mistake you are trying to point out
  
  vector<Particle> new_particles;

  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  // generate random starting index for resampling wheel
  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
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
  particle.associations= associations;
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