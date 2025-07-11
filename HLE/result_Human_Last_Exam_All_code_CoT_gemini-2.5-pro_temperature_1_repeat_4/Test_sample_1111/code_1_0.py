import random

def simulate(positions):
    """
    Simulates the particle activation process for a given set of initial positions.

    Args:
        positions: A sorted list of initial particle positions.

    Returns:
        The time T when the first particle hits the origin.
    """
    active_particles = [positions[0]]
    sleeping_positions = set(positions[1:])
    time = 0
    
    # We add a maximum time to prevent simulations for k<3 from running indefinitely,
    # as their expected time is infinite.
    MAX_TIME = 2000000 
    
    while True:
        time += 1
        if time > MAX_TIME:
            return MAX_TIME
            
        newly_activated_this_step = []
        next_positions = []

        for p in active_particles:
            # Each particle moves left or right with probability 1/2.
            move = 1 if random.random() < 0.5 else -1
            p_new = p + move
            
            # Check if a particle has reached the origin.
            if p_new == 0:
                return time
            
            next_positions.append(p_new)

        # After all particles have moved, check for new activations.
        for p_new in next_positions:
             if p_new in sleeping_positions:
                # A sleeping particle is activated.
                sleeping_positions.remove(p_new)
                newly_activated_this_step.append(p_new)

        # The list of active particles for the next step includes all previous
        # particles at their new locations, plus the newly activated ones.
        active_particles = next_positions + newly_activated_this_step

def main():
    """
    Runs simulations to compare the behavior of E[T] for k=2 and k=3.
    """
    N_TRIALS = 5000
    positions_k2 = [5, 10]
    positions_k3 = [5, 10, 15]
    
    # Checkpoints to print the running average
    checkpoints = sorted(list(set([100, 500, 1000, 2000, 5000, N_TRIALS])))

    # --- Simulation for k=2 ---
    print(f"--- Simulating for k=2 with initial positions {positions_k2} ---")
    print("We expect the average time to be unstable and grow, indicating infinite expectation.")
    total_time_k2 = 0
    for i in range(1, N_TRIALS + 1):
        total_time_k2 += simulate(positions_k2)
        if i in checkpoints:
            avg_time = total_time_k2 / i
            print(f"After {i} trials, average T for k=2: {avg_time:.2f}")

    print("\n" + "="*50 + "\n")

    # --- Simulation for k=3 ---
    print(f"--- Simulating for k=3 with initial positions {positions_k3} ---")
    print("We expect the average time to converge to a finite value.")
    total_time_k3 = 0
    for i in range(1, N_TRIALS + 1):
        total_time_k3 += simulate(positions_k3)
        if i in checkpoints:
            avg_time = total_time_k3 / i
            print(f"After {i} trials, average T for k=3: {avg_time:.2f}")

if __name__ == '__main__':
    main()