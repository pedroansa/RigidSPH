#include "simulation.hpp"


using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Fill it with laplacian of W_viscosity
    vec3 r = p_i - p_j;
    float r_norm = norm(r);
    float par1 = 45.0f / (3.14159f * std::pow(h, 6));
    float par2 = (h - r_norm);

    assert_cgp_no_msg(r_norm <= h);

    return par1 * par2;
    
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    vec3 r = p_i - p_j;
    float r_norm = norm(r);

    float par1 = -45.0f / (3.14159f * std::pow(h, 6));
    float par2 = (h - r_norm) * (h - r_norm);

    assert_cgp_no_msg(r_norm <= h);
    return par1 * par2 * r/r_norm;
   
    
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    //std::cout << r << std::endl;
    assert_cgp_no_msg(r <= h);
    return 315.0 / (64.0 * 3.14159f * std::pow(h, 9)) * std::pow(h * h - r * r, 3.0f);


}


void update_density(numarray<particle_element>& particles, float h, float m)
{
    // Compute the density value (particles[i].rho) at each particle position
    //  rho_i = \sum_j m W_density(pi,pj)
    int const N = particles.size();
    for (int i = 0; i < N; ++i) {
        float rho_i = 0.0f; // 
        particle_element p1 = particles[i];

        // Iterate over all particles to compute the density contribution from each neighbor
        for (int j = 0; j < N; ++j) {

            particle_element p2 = particles[j];
            float const r = norm(p1.p - p2.p);
            //std::cout << W << std::endl;
            if (r < h) {
                float W = W_density(p1.p, p2.p, h);
                rho_i += m * W; 
            }

        }

        particles[i].rho = rho_i;
    }


}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
void update_force(numarray<particle_element>& particles, float h, float m, float nu)
{
    // gravity
    const int N = particles.size();
    for (int i = 0; i < N; ++i)
        particles[i].f = m * vec3{ 0,-9.81f,0 };

    // For all particles i
    //   Compute F_pressure
    //   Compute F_viscosity
    for (int i = 0; i < N; ++i) {
        vec3 f_pressure = vec3{ 0, 0, 0 };
        vec3 f_viscosity = vec3{ 0, 0, 0 };
        particle_element& p1 = particles[i];

        // Calcular for�a de press�o
        vec3 aux = vec3(0.0f,0.0f,0.0f);
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                particle_element& p2 = particles[j];
                if (norm(p1.p - p2.p) <= h) {
                    vec3 W = W_gradient_pressure(p1.p, p2.p, h);

                    aux += m * 0.5f*((p1.pressure + p2.pressure) / p2.rho) * W;
                }
            }
        }
        f_pressure = -(m / p1.rho) * aux;

        // Viscocity
        aux = vec3(0.0f, 0.0f, 0.0f);
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                particle_element& p2 = particles[j];
                if (norm(p1.p - p2.p) <= h) {
                    float w = W_laplacian_viscosity(p1.p, p2.p, h);

                    aux += m * ((p2.v - p1.v) / p2.rho) * w;
                }
            }
        }
        f_viscosity = m * nu * aux;

        particles[i].f += (f_pressure + f_viscosity);
    }
}

/*void apply_force_at_point(rigid_body& square, const cgp::vec3& force, const cgp::vec3& point) {
    // Calculate torque by taking the cross product of the force and the vector from the center of mass to the point
    cgp::vec3 r = point - square.position;
    cgp::vec3 torque = cross(r, force);

    // Apply linear force
    square.force += force;

    // Apply torque
    square.torque += torque;
}*/

/*float cubic_spline_kernel(float r, float h) {
    const float coeff = 1.0f / (7.0f * 3.14159f * h * h * h);
    if (r <= h) {
        float q = 1 - r / h;
        return coeff * q * q * q;
    }
    return 0.0f;
}

// Calculate volume of a particle based on its neighbors
float calculate_particle_volume(const vec3& particle_position, const numarray<particle_element>& particles, float h) {
    float volume = 0.0f;
    const int N = particles.size();
    for (int i = 0; i < N; ++i) {
        const vec3& neighbor_position = particles[i].p;
        float distance = norm(particle_position - neighbor_position);
        volume += cubic_spline_kernel(distance, h);
    }
    return volume;
}

void update_buoyancy(numarray<particle_element>& particles, rigid_body& square, float fluid_density, sph_parameters_structure const& sph_parameters)
{
    // Calculate the displaced fluid volume
    float displaced_volume = 0.0f;
    float half_width = square.square_size.x / 2.0f;
    float half_height = square.square_size.y / 2.0f;

    // Iterate over all particles to check if they're inside the square
    for (const auto& particle : particles) {
        vec3 particle_position = particle.p;
        // Check if the particle is inside the bounding box of the square
        if (particle_position.x >= square.position.x - half_width &&
            particle_position.x <= square.position.x + half_width &&
            particle_position.y >= square.position.y - half_height &&
            particle_position.y <= square.position.y + half_height) {
            // Particle is inside the square, add its volume to the displaced volume
            displaced_volume += calculate_particle_volume(particle_position, particles, sph_parameters.h)/20000.0f; // Call calculate_particle_volume // You need to define how to calculate particle volume
        }
    }

    // Calculate the buoyancy force
    float buoyant_force_magnitude = fluid_density * 9.81f * displaced_volume;

    // Apply the buoyant force in the opposite direction of gravity
    cgp::vec3 buoyant_force = { 0, buoyant_force_magnitude, 0 };

    // Apply buoyant force to the square
    square.force += buoyant_force;
}*/

void update_force_solid(rigid_body& square) {
    const float gravitational_constant = 9.81f; // m/s^2 (adjust as needed)

    // Reset force to zero
    square.force = { 0.0f, 0.0f, 0.0f };
    //square.torque = { 0.0f, 0.0f, 0.0f };

    // Calculate gravitational force (F = m * g)
    cgp::vec3 gravitational_force = { 0, -square.mass * gravitational_constant, 0 };

    // Apply gravitational force
    square.force += gravitational_force;

    /*cgp::vec3 force_at_corner = {10.0f, 0.0f, 0.0f};
    cgp::vec3 corner_position = square.position + cgp::vec3(square.square_size.x / 2, square.square_size.y / 2, 0); // Adjust for the actual corner position
    //apply_force_at_point(square, force_at_corner, corner_position);*/
}

void simulate(float dt, numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters, rigid_body& square)
{

	// Update values
    update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces

    update_force_solid(square);

    float submerged_volume = 0.0f;
    const float particle_radius = sph_parameters.h; // Assuming particles have a spherical interaction radius
    const int N = particles.size();
    for (int i = 0; i < N; ++i) {
        const vec3& p = particles[i].p;
        // Check if particle is within the bounds of the square
        if (p.x >= square.position.x - square.square_size.x / 2 &&
            p.x <= square.position.x + square.square_size.x / 2 &&
            p.y >= square.position.y - square.square_size.y / 2 &&
            p.y <= square.position.y + square.square_size.y / 2) {
            // Compute submerged volume using particle-based approximation
            submerged_volume += 100.0f / 3.0f * 3.14159f * std::pow(particle_radius, 3);
        }
    }

    // Buoyancy force = displaced volume * density of fluid * gravity
    const float density_fluid = sph_parameters.rho0; // Assuming fluid density equals reference density
    const float buoyancy_force = submerged_volume * density_fluid * 9.81f;

    // Apply buoyancy force
    square.force += vec3(0, buoyancy_force, 0);
    std::cout << square.force << std::endl;
    // Update solid square position and velocity
    square.velocity += (1 - 0.005f) * dt * (square.force / square.mass);
    square.position += dt * square.velocity;

    /*const float fluid_density = 100.0f;
    update_buoyancy(particles, square, fluid_density, sph_parameters);

    square.velocity = square.velocity + dt * (square.force / square.mass); // Update velocity
    square.position = square.position + dt * square.velocity; // Update position*/

    //----------------------------------------

    // Update angular motion
    /*square.angular_acceleration = square.torque / square.mass; // Assuming constant mass for simplicity
    square.angular_velocity += dt * square.angular_acceleration;

    // Update mesh orientation based on angular velocity
    float angle = norm(square.angular_velocity) * dt; // Calculate the angle of rotation
    cgp::vec3 axis = normalize(square.angular_velocity); // Calculate the axis of rotation

    // Construct the quaternion representing the rotation
    cgp::quaternion rotation_quaternion = cgp::quat_from_axis_angle(axis, angle);

    // Update the rotation of the mesh
    square.mesh.model.rotation = rotation_quaternion;*/

    // Update mesh translation to match the position of the square
    square.mesh.model.translation = square.position;

    float half_width = square.square_size.x / 4.0f;
    float half_height = square.square_size.y / 4.0f;

    // Collision handling for the square
    float const epsilon = 1e-3f;
    if (square.position.y - half_height < -1) {
        square.position.y = -1 + half_height + epsilon * rand_uniform();
        square.velocity.y *= -0.5f;
    }

    if (square.position.x - half_width < -1) {
        square.position.x = -1 + half_width + epsilon * rand_uniform();
        square.velocity.x *= -0.5f;
    }

    if (square.position.x + half_width > 1) {
        square.position.x = 1 - half_width - epsilon * rand_uniform();
        square.velocity.x *= -0.5f;
    }

	// Numerical integration
	float const damping = 0.005f;
    //int const N = particles.size();
	float const m = sph_parameters.m;
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}

	// Collision
    for(int k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
    }

}