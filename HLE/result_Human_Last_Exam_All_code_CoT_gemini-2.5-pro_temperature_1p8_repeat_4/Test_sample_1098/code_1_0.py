import random

def simulate_brownian_ratchet(temperature, total_particles=1000, max_steps=10000):
    """
    Simulates a simplified Brownian ratchet to show the role of temperature.

    Args:
        temperature (float): The temperature of the system. A value of 0 implies no thermal motion.
        total_particles (int): The total number of gas particles.
        max_steps (int): The maximum number of simulation steps to run.
    """
    print(f"--- Running Simulation ---")
    print(f"Parameter Tested: Temperature = {temperature} K")

    # Start with an equal number of particles in each compartment
    compartment_a = total_particles // 2
    compartment_b = total_particles - compartment_a

    print(f"Initial State -> Compartment A: {compartment_a}, Compartment B: {compartment_b}")

    # The probability of a particle crossing the one-way door is
    # directly related to its thermal motion, which is a function of temperature.
    # If T=0, there's no motion, so the probability of moving is 0.
    # If T>0, there's motion, so there's a non-zero probability.
    move_probability = 0.01 if temperature > 0 else 0

    if move_probability == 0:
        print("At temperature = 0, there is no thermal motion. Particles will not move.")
    else:
        print("At temperature > 0, particles are in random motion.")

    # Run the simulation for a number of steps
    for i in range(max_steps):
        if compartment_a == 0:
            print(f"\nProcess completed at step {i}.")
            break

        # In each step, a fraction of particles in compartment A, determined
        # by the move_probability, will move to B due to random motion.
        # This is a simplification of a binomial process.
        moved_count = 0
        for _ in range(compartment_a):
            if random.random() < move_probability:
                moved_count += 1
        
        if moved_count > 0:
            compartment_a -= moved_count
            compartment_b += moved_count

    print(f"Final State -> Compartment A: {compartment_a}, Compartment B: {compartment_b}")
    print(f"--------------------------\n")

# --- Main Execution ---

# Case 1: Test with Temperature at absolute zero.
# At T=0, there is no kinetic energy, no motion, and no particles can cross the door.
simulate_brownian_ratchet(temperature=0)

# Case 2: Test with a positive Temperature (e.g., room temperature ~298 K).
# At T>0, thermal motion causes particles to cross from A to B until all are trapped.
simulate_brownian_ratchet(temperature=298)