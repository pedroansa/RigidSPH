#pragma once

#include "cgp/cgp.hpp"

struct rigid_body
{
	cgp::mesh_drawable mesh; // Visual representation
	cgp::vec3 position;      // Position of the rigid body's center
	cgp::vec3 velocity;      // Linear velocity
	cgp::vec3 acceleration;  // Linear acceleration
	cgp::vec3 square_size;

	cgp::vec3 force; // Force
	cgp::vec3 torque;// Torque

	cgp::vec3 angular_velocity;        // Angular velocity
	cgp::vec3 angular_acceleration;    // Angular acceleration

	float mass;

	// Constructor
	rigid_body() : position(0, 0, 0), velocity(0, 0, 0), acceleration(0, 0, 0), mass(1.0f),
		angular_velocity(0, 0, 0), angular_acceleration(0, 0, 0), torque(0, 0, 0) // Initialize angular velocity, angular acceleration, and torque to zero
	{
		square_size = { 0.2f, 0.2f, 0.0f };

	}

	// Function to update the physics of the rigid body
	/*void update(float dt)
	{
		acceleration.y -= 9.8f;
		velocity += acceleration * dt; // Update velocity with acceleration
		position += velocity * dt;     // Update position with velocity

		// Update the mesh position to match the rigid body's position
		mesh.model.translation = position;

		// Reset acceleration if it's being applied anew each frame
		acceleration = cgp::vec3(0, 0, 0);
	}*/
};

// SPH Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},rho(0),pressure(0) {}
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.12f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

     // Total mass of a particle (consider rho0 h^2)
    float m = rho0*h*h;

    // viscosity parameter
    float nu = 0.02f;   
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
};


void simulate(float dt, cgp::numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters, rigid_body& square);