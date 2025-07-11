import numpy as np
import random

def run_hitting_probability_simulation():
    """
    Simulates a 3D random walk to estimate the hitting probability of a transient set.
    """
    d = 3
    # A single point is a classic example of a transient set in d>=3.
    # We represent it as a set of tuples for efficient lookup.
    A_set = {(1, 0, 0)}
    a_point = np.array(list(A_set)[0])

    # Simulation parameters
    num_trials = 20000  # Number of random walks to simulate for each starting point
    max_steps = 50000  # Max steps for a single walk before assuming it escaped

    # Define the possible moves (steps) for a simple random walk on Z^d
    moves = []
    for i in range(d):
        for sign in [-1, 1]:
            move = np.zeros(d, dtype=int)
            move[i] = sign
            moves.append(move)

    # Start points with increasing distance from the set A
    start_points = [np.array([k, 0, 0]) for k in range(2, 12)]

    print("--- Simulating Hitting Probability for a Transient Set A = {(1, 0, 0)} in 3D ---")
    print(f"This demonstrates that for a transient set, P(hit A) -> 0 as distance increases.\n")
    
    # Store results to show theoretical relationship
    results = []

    for start_pos in start_points:
        hit_count = 0
        for _ in range(num_trials):
            current_pos = np.copy(start_pos)
            for _ in range(max_steps):
                # Take a random step
                step = random.choice(moves)
                current_pos += step

                # Check if the walker has hit the set A
                if tuple(current_pos) in A_set:
                    hit_count += 1
                    break
        
        hit_prob = hit_count / num_trials
        dist = np.linalg.norm(start_pos - a_point)
        results.append((dist, hit_prob))

        print(f"For starting point x = {tuple(start_pos)}, distance from A = {dist:.1f}:")
        # "output each number in the final equation!"
        print(f"  P(hit A) = {hit_count} / {num_trials} = {hit_prob:.4f}")

    print("\n--- Analysis ---")
    print("Theory predicts that P(hit A) is proportional to 1/distance for a point set in 3D.")
    print("We can check this by observing the product of probability and distance:")
    
    # Calculate Product = P(hit) * distance
    for dist, prob in results:
        product = prob * dist
        print(f"  Dist: {dist:4.1f}, P(hit): {prob:.4f}  => P(hit) * Dist = {product:.4f}")

    print("\nThe simulation confirms that the hitting probability decreases with distance and is never 1,")
    print("which is the behavior expected of a transient set, and opposite to the condition in the problem.")

if __name__ == '__main__':
    run_hitting_probability_simulation()
<<<No, such a set A cannot be transient.>>>