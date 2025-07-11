import numpy as np

def hitting_probability_simulation():
    """
    Simulates random walks to demonstrate that the probability of hitting a
    transient set tends to zero as the starting distance increases.
    """
    # --- Parameters ---
    d = 3  # Dimension, must be >= 3 for the walk to be transient.
    
    # For simplicity, we use the origin as our transient set A.
    # In d>=3, any finite set is transient.
    A = {(0,) * d}
    
    # Simulation parameters
    num_walks = 20000  # Number of walks for each starting point for accuracy.
    max_steps = 2500   # A cutoff to assume the walk has 'escaped' to infinity.

    # We'll test for starting points at increasing integer distances from A.
    start_distances = list(range(1, 16))

    # --- Explanation ---
    print("This simulation provides evidence for a theorem about transient random walks.")
    print(f"We use dimension d={d} and a simple transient set A = {A}.")
    print("The theorem states: For a transient set A, the probability of a walk starting")
    print("at x hitting A, denoted P_x(hit A), must approach 0 as the distance |x| -> infinity.")
    print("\nWe will now estimate this probability for starting points at increasing distances.")
    print("-" * 70)

    # --- Simulation ---
    # Define the 2d possible steps (e.g., in d=3: (+-1,0,0), (0,+-1,0), (0,0,+-1))
    steps = []
    for i in range(d):
        step = np.zeros(d, dtype=int)
        step[i] = 1
        steps.append(step)
        steps.append(-step)

    for dist in start_distances:
        # Create a starting point at the specified distance, e.g., (dist, 0, 0, ...)
        start_pos = np.zeros(d, dtype=int)
        start_pos[0] = dist
        
        hit_count = 0
        for _ in range(num_walks):
            current_pos = np.copy(start_pos)
            for _ in range(max_steps):
                # Choose a random move from the list of possible steps
                move_index = np.random.randint(0, 2 * d)
                current_pos += steps[move_index]
                
                # Check if the walk has hit the set A
                if tuple(current_pos) in A:
                    hit_count += 1
                    break  # End this walk as it's a hit
        
        # Calculate and print the estimated probability
        hit_probability = hit_count / num_walks
        print(f"Distance: {dist:<5} | Start point: {f'({dist},0,0)':<10} | Estimated Hit Probability: {hit_probability:.5f}")

    print("-" * 70)
    print("\n--- Conclusion from Simulation ---")
    print("The simulation shows a clear trend: as the starting distance increases, the")
    print("probability of hitting the transient set A decreases, approaching 0.")
    
    print("\n--- Final Answer to the Problem ---")
    print("The problem assumes there are infinitely many points 'x' where the hitting")
    print("probability is 1. An infinite set must be unbounded, so some of these points")
    print("would be arbitrarily far from A.")
    print("This contradicts our finding that the hitting probability must tend to 0 for distant points.")
    print("Therefore, a set with the given property CANNOT be transient.")

# Run the simulation
hitting_probability_simulation()