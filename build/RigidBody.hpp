#pragma once
#include "cgp/cgp.hpp"

class RigidBody
{
public:
	cgp::mesh_drawable mesh; // Visual representation
	cgp::vec3 position;      // Position of the rigid body's center
	cgp::vec3 velocity;      // Linear velocity
	cgp::vec3 acceleration;  // Linear acceleration

	// Constructor
	RigidBody() : position(0, 0, 0), velocity(0, 0, 0), acceleration(0, 0, 0)
	{



	}

	// Function to update the physics of the rigid body
	void update(float dt)
	{
		acceleration.y -= 9.8f;
		velocity += acceleration * dt; // Update velocity with acceleration
		position += velocity * dt;     // Update position with velocity

		// Update the mesh position to match the rigid body's position
		mesh.model.translation = position;

		// Reset acceleration if it's being applied anew each frame
		acceleration = cgp::vec3(0, 0, 0);
	}
};

