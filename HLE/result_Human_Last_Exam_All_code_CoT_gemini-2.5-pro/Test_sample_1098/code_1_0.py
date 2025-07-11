import random

def simulate_maxwell_demon_door():
    """
    Simulates an idealized Maxwell's demon apparatus with a perfect one-way door.
    """
    # Initial conditions: Equal amounts of gas in each compartment.
    # We represent the "amount" by the number of particles.
    particles_left = 50
    particles_right = 50
    total_particles = particles_left + particles_right

    print("This simulation demonstrates the behavior of a 'perfect' one-way door.")
    print("In this ideal case, particles can only move from the left to the right compartment.")
    print("-" * 40)
    print(f"Initial State | Left: {particles_left}, Right: {particles_right}")
    print("-" * 40)

    time_step = 0
    # The simulation continues until the left compartment is empty.
    while particles_left > 0:
        time_step += 1

        # In any given moment, a particle from either side might hit the door.
        # The chance of a particle hitting the door is proportional to the number
        # of particles on that side. We simulate this by checking if a random
        # particle chosen from the total is on the left side.
        if random.uniform(0, total_particles) < particles_left:
            # A particle from the left hits the door and successfully passes through.
            particles_left -= 1
            particles_right += 1

        # If a particle from the right hits the door, it is blocked.
        # So, we do not need to model this event as it results in no change.

        # Print the state periodically to observe the progression.
        if time_step % 50 == 0:
            print(f"Time Step {time_step:4} | Left: {particles_left:2}, Right: {particles_right:2}")

    print("-" * 40)
    print(f"Final State   | Left: {particles_left}, Right: {particles_right} at Time Step {time_step}")
    print("-" * 40)
    print("\nExplanation of the physical requirement:")
    print("The simulation shows all gas moving to one side. This relied on a 'perfect' door.")
    print("A real door at the same temperature as the gas would jiggle due to its own thermal energy,")
    print("allowing particles to leak backward, preventing this outcome.")
    print("To make the door truly one-way, it must be at a different (much lower) temperature than the gas.")
    print("Therefore, TEMPERATURE is the critical experimental parameter.")

simulate_maxwell_demon_door()