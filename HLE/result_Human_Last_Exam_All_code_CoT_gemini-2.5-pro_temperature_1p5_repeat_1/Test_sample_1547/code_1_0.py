import random
import math

def estimate_hitting_probability():
    """
    Estimates the hitting probability of a transient set A for a simple
    random walk in d=3, starting from points at increasing distances.
    """
    d = 3
    # A is the origin, a transient set for d>=3.
    A = {(0, 0, 0)}

    # Number of simulations for each starting point.
    num_simulations = 20000

    print("Simulating SRW in 3D to estimate hitting probability of the origin.")
    print("------------------------------------------------------------------")
    print("Dist | Start Point   | Estimated Hitting Probability")
    print("------------------------------------------------------------------")

    # List of basis vectors for random steps.
    basis_vectors = []
    for i in range(d):
        pos_vec = [0] * d
        pos_vec[i] = 1
        basis_vectors.append(tuple(pos_vec))
        neg_vec = [0] * d
        neg_vec[i] = -1
        basis_vectors.append(tuple(neg_vec))

    # Test for starting points at various distances from the origin.
    for dist in range(1, 11):
        # Start at (dist, 0, 0).
        start_pos = tuple([dist] + [0] * (d - 1))

        # Max steps should scale with distance squared.
        max_steps = 5 * dist**2 + 200
        
        hit_count = 0
        for _ in range(num_simulations):
            # Start a new walk.
            current_pos = list(start_pos)
            
            for _ in range(max_steps):
                # Take one random step.
                step = random.choice(basis_vectors)
                for i in range(d):
                    current_pos[i] += step[i]

                # Check if the walk hits the set A.
                if tuple(current_pos) in A:
                    hit_count += 1
                    break
        
        # Calculate and print the estimated probability.
        probability = hit_count / num_simulations
        print(f"{dist:<4} | {str(start_pos):<13} | {probability:.4f}")

    print("------------------------------------------------------------------")
    print("\nAs observed, the probability of hitting the origin decreases")
    print("as the starting distance increases, approaching 0.")
    print("This illustrates that it's impossible for the probability to be 1")
    print("for infinitely many (and thus arbitrarily far) points.")

if __name__ == '__main__':
    estimate_hitting_probability()