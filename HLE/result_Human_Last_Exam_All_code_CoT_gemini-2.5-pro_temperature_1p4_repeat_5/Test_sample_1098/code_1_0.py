import random

def run_simulation():
    """
    Simulates a Maxwell's demon setup to demonstrate the importance of temperature.
    """
    total_particles = 1000
    simulation_steps = 5000

    print("This simulation demonstrates why a temperature difference is required to trap gas on one side of a compartment.")
    print("The 'one-way door' is a Brownian ratchet, whose function depends on its temperature relative to the gas.\n")

    # --- Scenario 1: System at Uniform Temperature ---
    print("--- SCENARIO 1: Uniform Temperature (Door and Gas are the same T) ---")
    print("The door has thermal energy, causing it to jiggle randomly. This allows particles to pass in BOTH directions.")
    
    # Initial state
    particles_A = total_particles // 2
    particles_B = total_particles - particles_A

    for _ in range(simulation_steps):
        # A particle has a chance to move from A to B
        if random.random() < particles_A / total_particles:
            particles_A -= 1
            particles_B += 1
        
        # A particle also has a chance to move from B to A because the door jiggles back
        if random.random() < particles_B / total_particles:
            particles_B -= 1
            particles_A += 1

    print("Final State at Uniform Temperature:")
    print(f"Compartment A has {particles_A} particles.")
    print(f"Compartment B has {particles_B} particles.")
    # The "final equation" represents the final distribution
    print(f"Final Equation: {particles_A} (A) + {particles_B} (B) = {particles_A + particles_B} (Total)")
    print("Result: The system remains near equilibrium, with gas distributed between both sides.\n")


    # --- Scenario 2: Temperature Difference (Door is colder than Gas) ---
    print("--- SCENARIO 2: Temperature Difference (Door << Gas T) ---")
    print("The door is very cold and has no thermal energy to jiggle back. It acts as a perfect one-way gate from A to B.")
    
    # Reset to initial state
    particles_A = total_particles // 2
    particles_B = total_particles - particles_A

    # In this simulation, we'll model movement from A to B only
    # We increase the number of steps to ensure completion
    for _ in range(simulation_steps * 2):
        if particles_A > 0:
            # We model a particle hitting the door from side A and passing through
            # The probability is higher to speed up the simulation of the concept
            if random.random() < 0.75:
                particles_A -= 1
                particles_B += 1
        # No movement from B to A is allowed

    print("Final State with a Temperature Difference:")
    print(f"Compartment A has {particles_A} particles.")
    print(f"Compartment B has {particles_B} particles.")
    # The "final equation" represents the final distribution
    print(f"Final Equation: {particles_A} (A) + {particles_B} (B) = {particles_A + particles_B} (Total)")
    print("Result: All gas is trapped in one compartment, which is only possible with a temperature difference.\n")

    print("Conclusion: The essential experimental parameter is TEMPERATURE.")

if __name__ == '__main__':
    run_simulation()