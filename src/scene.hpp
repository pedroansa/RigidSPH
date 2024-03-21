#pragma once


#include "cgp/cgp.hpp"
#include "environment.hpp"

#include "simulation/simulation.hpp"

using cgp::mesh_drawable;


struct gui_parameters {
	bool display_color = true;
	bool display_particles = true;
	bool display_radius = false;
};

struct rigid_body
{
	cgp::mesh_drawable mesh; // Visual representation
	cgp::vec3 position;      // Position of the rigid body's center
	cgp::vec3 velocity;      // Linear velocity
	cgp::vec3 acceleration;  // Linear acceleration
	vec3 square_size;

	// Constructor
	rigid_body() : position(0, 0, 0), velocity(0, 0, 0), acceleration(0, 0, 0)
	{
		square_size = { 0.2f, 0.2f, 0.0f };

		
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


// The structure of the custom scene
struct scene_structure : cgp::scene_inputs_generic {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	camera_controller_orbit camera_control;
	camera_projection_orthographic camera_projection;
	window_structure window;

	mesh_drawable global_frame;          // The standard global frame
	environment_structure environment;   // Standard environment controler
	input_devices inputs;                // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                  // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	sph_parameters_structure sph_parameters; // Physical parameter related to SPH
	cgp::numarray<particle_element> particles;      // Storage of the particles
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle
	cgp::curve_drawable curve_visual;   // Circle used to display the radius h of influence

	cgp::grid_2D<cgp::vec3> field;      // grid used to represent the volume of the fluid under the particles
	cgp::mesh_drawable field_quad; // quad used to display this field color

	rigid_body rigid_body_square;


	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();    // Standard initialization to be called before the animation loop
	void display_frame(); // The frame display to be called within the animation loop
	void display_gui();   // The display of the GUI, also called within the animation loop

	void initialize_sph();

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

};





