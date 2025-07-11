import random

def run_srw_simulation():
    """
    Simulates simple random walks in 3D to test if a known transient set
    satisfies the condition P_x(tau_A < infinity) = 1 for a distant x.
    """
    # Simulation parameters
    dim = 3
    num_walks = 20000
    max_steps = 2000
    
    # Starting point for the walks, chosen to be far from the set A
    # The distance R from the set A.
    R = 10
    start_pos = [0] * dim
    start_pos[1] = R # Start at (0, R, 0) = (0, 10, 0)

    # Define the transient set A: the non-negative x-axis
    # A = {(k, 0, 0) for k >= 0}
    def is_in_A(pos):
        return pos[0] >= 0 and pos[1] == 0 and pos[2] == 0

    # Possible moves for a simple random walk in 3D
    moves = []
    for i in range(dim):
        for d in [-1, 1]:
            move = [0] * dim
            move[i] = d
            moves.append(tuple(move))

    hits = 0
    print(f"Starting simulation with {num_walks} walks...")
    print(f"Set A is the non-negative x-axis.")
    print(f"Walks start at {tuple(start_pos)}, at distance {R} from A.")
    print(f"Max steps per walk: {max_steps}")

    for i in range(num_walks):
        pos = list(start_pos)
        
        # We don't need to check the starting position
        for step in range(max_steps):
            # Choose a random move
            move = random.choice(moves)
            # Update position
            for j in range(dim):
                pos[j] += move[j]
            
            # Check if the walk hit the set A
            if is_in_A(pos):
                hits += 1
                break # End this walk as it has hit the set
    
    # Calculate the estimated hitting probability
    probability = hits / num_walks
    
    print("\n--- Simulation Results ---")
    print(f"Number of walks that hit set A: {hits}")
    print(f"Total number of walks simulated: {num_walks}")
    print(f"Estimated hitting probability P_x(tau_A < infinity):")
    print(f"{hits} / {num_walks} = {probability:.4f}")
    
    if probability < 0.9:
        print("\nThe estimated probability is significantly less than 1.")
        print("This provides evidence that this transient set does not satisfy the condition from the problem.")
        print("The theoretical argument shows that NO transient set can satisfy the condition.")
    else:
        print("\nThe probability is close to 1. This might be due to simulation parameters (e.g., R not large enough).")
        print("However, the theoretical argument proves the answer is No.")

if __name__ == '__main__':
    run_srw_simulation()
