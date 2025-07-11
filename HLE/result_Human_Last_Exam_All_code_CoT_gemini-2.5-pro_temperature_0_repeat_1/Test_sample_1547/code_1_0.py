import random
import math

def run_simulation():
    """
    Main function to run and print the simulation results.
    """

    # A simple random walk in 3D takes a step in one of 6 directions.
    DIRECTIONS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

    def single_walk(start_pos, is_in_set_A, max_steps):
        """
        Performs a single random walk in 3D starting from start_pos.
        The walk is for n>=1, so we check for a hit after the first step.
        Returns True if the walk hits set A within max_steps, False otherwise.
        """
        pos = list(start_pos)
        for _ in range(max_steps):
            step = random.choice(DIRECTIONS)
            pos[0] += step[0]
            pos[1] += step[1]
            pos[2] += step[2]
            if is_in_set_A(tuple(pos)):
                return True
        return False

    def estimate_hitting_prob(start_pos, is_in_set_A, num_walks, max_steps):
        """
        Estimates the probability of a random walk hitting set A by running
        many simulations.
        """
        if is_in_set_A(start_pos):
            return 1.0
            
        hits = 0
        for _ in range(num_walks):
            if single_walk(start_pos, is_in_set_A, max_steps):
                hits += 1
        return hits / num_walks

    # --- Define the sets for our simulation ---

    # 1. A transient set: a single point, the origin. Any finite set is transient.
    def is_in_A_transient(pos):
        return pos == (0, 0, 0)

    # 2. A recurrent set: the hyperplane z=0. This set satisfies the problem's
    #    condition, as P_x(tau_A < inf) = 1 for all x.
    def is_in_A_recurrent(pos):
        return pos[2] == 0

    # --- Simulation Parameters ---
    NUM_WALKS = 10000  # Number of simulations for each data point
    BASE_MAX_STEPS = 200 # Base for maximum steps in a single walk

    print("--- Illustration for the Theoretical Argument ---")
    print("\nPart 1: Simulating a TRANSIENT set A = {(0,0,0)}")
    print("="*60)
    print("Theory: The hitting probability P_x(tau_A < inf) should approach 0 as |x| increases.")
    print("We estimate this probability for starting points x at increasing distances.")
    print("-" * 60)
    for k in range(1, 8):
        start_point = (k, 0, 0)
        # The time to hit can be long, so we increase max_steps for farther points.
        # The factor k*k is related to the diffusive scaling.
        max_s = BASE_MAX_STEPS * k * k
        prob = estimate_hitting_prob(start_point, is_in_A_transient, NUM_WALKS, max_s)
        print(f"Start x = {start_point}, |x|={k}: Estimated hitting probability = {prob:.4f}")

    print("\n\nPart 2: Simulating a RECURRENT set A = {plane z=0}")
    print("="*60)
    print("Theory: The hitting probability P_x(tau_A < inf) is 1 for all x.")
    print("Our simulation should yield values close to 1, limited only by max_steps.")
    print("-" * 60)
    for k in range(1, 8):
        start_point = (0, 0, k)
        max_s = BASE_MAX_STEPS * k * k
        prob = estimate_hitting_prob(start_point, is_in_A_recurrent, NUM_WALKS, max_s)
        print(f"Start x = {start_point}, dist={k}: Estimated hitting probability = {prob:.4f}")
        
    print("\n--- End of Illustration ---")
    print("The simulation supports the theory: The property P_x(tau_A < inf) = 1 for infinitely many x")
    print("is characteristic of recurrent sets, not transient ones.")

run_simulation()