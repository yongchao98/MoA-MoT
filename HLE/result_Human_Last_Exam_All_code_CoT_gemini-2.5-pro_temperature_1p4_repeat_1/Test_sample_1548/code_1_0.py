import numpy as np

def solve_task():
    """
    This script addresses the question about transient sets for the h-transformed
    simple random walk on Z^2. It provides both a numerical calculation
    to build intuition and a theoretical explanation for the final answer.
    """

    # --- Part 1: Numerical Calculation ---
    # We analyze the transition probabilities from a point x = (k, 0) on the x-axis,
    # for a large k, to illustrate the behavior of the walk.

    # The potential kernel a(x) for SRW on Z^2 is approximated by:
    # a(x) approx (2/pi) * log(||x||) + C
    # where C is a constant (related to Euler's constant gamma and log(8)).
    GAMMA = 0.5772156649
    LOG_8 = np.log(8)
    C = (2 * GAMMA + LOG_8) / np.pi

    def potential_kernel_approx(point):
        """Approximation of the potential kernel a(z) for a point z."""
        x, y = point
        # The state space is Z^2 - {0}. The potential kernel a(0) is 0.
        if x == 0 and y == 0:
            return 0.0
        
        norm = np.sqrt(x**2 + y**2)
        return (2 / np.pi) * np.log(norm) + C

    k = 100
    start_point = (k, 0)
    print(f"--- Numerical Analysis ---")
    print(f"Let's analyze the walk's transitions from the point x = {start_point}.")

    # The four neighbors of (k, 0)
    neighbors = {
        "right": (k + 1, 0),
        "left": (k - 1, 0),
        "up": (k, 1),
        "down": (k, -1)
    }
    
    # The value of the harmonic function h at the starting point x
    h_x = potential_kernel_approx(start_point)

    print(f"The value of the harmonic function h(x) at x is: a({start_point}) ~ {h_x:.4f}")
    print("\nCalculating transition probabilities P_hat(x,y) = (1/4) * h(y)/h(x):")
    print("-" * 55)
    print(f"{'Neighbor (y)':<15} | {'h(y)':<10} | {'Probability P_hat(x,y)':<25}")
    print("-" * 55)
    
    total_prob = 0
    
    for name, pos in neighbors.items():
        h_y = potential_kernel_approx(pos)
        # For the h-transform, the new transition probability is P_hat(x,y) = P(x,y) * h(y)/h(x).
        # For SRW on Z^2, the standard transition prob P(x,y) is 1/4 for all neighbors.
        prob = (1 / 4) * h_y / h_x
        total_prob += prob
        print(f"{name:<15} | {h_y:<10.4f} | {prob:.8f}")
    
    print("-" * 55)
    print(f"Sum of calculated probabilities = {total_prob:.8f}")
    print("(Note: The sum is very close to 1. The small deviation is because our\n"
          " h(x) is an approximation. For the true potential kernel, which is perfectly\n"
          " harmonic, this sum would be exactly 1.)")
    print("\n")
    
    # --- Part 2: Theoretical Explanation ---
    print("--- Theoretical Answer ---")
    print("The question is: For the h-transformed SRW on Z^2, is it true that every transient set must necessarily be finite?")
    print("\nThe answer is NO.\n")
    
    print("Justification:")
    print("1. The h-transform with the potential kernel a(x) conditions the random walk to escape to infinity.")
    print("   A key result from the theory of random walks states that the direction of this process converges.")
    print("   Specifically, if S_hat_n is the position of the walk at step n, the vector S_hat_n / ||S_hat_n||")
    print("   converges almost surely to a random vector V on the unit circle.")
    
    print("\n2. Because the harmonic function a(x) is (asymptotically) radially symmetric, the limiting")
    print("   direction V is uniformly distributed on the circle. This means the walk's path will, with")
    print("   probability 1, eventually follow a straight line in a randomly chosen direction.")
          
    print("\n3. Let's test the claim with a counterexample. Consider the set A equal to the positive x-axis:")
    print("   A = { (k, 0) | k is an integer and k >= 1 }")
    print("   This set A is clearly an infinite set.")
          
    print("\n4. For a path to visit the set A infinitely many times, it must keep returning to the x-axis.")
    print("   This would imply that the limiting direction of the path must be along the x-axis (i.e., an angle of 0 or pi).")
    
    print("\n5. However, the probability that the limiting angle takes any specific value (like 0 or pi) is zero,")
    print("   because the distribution of the angle is continuous (uniform).")
    
    print("\n6. Therefore, with probability 1, the walk will only visit the set A a finite number of times before")
    print("   departing towards its final direction forever. This means the set A is a transient set.")
          
    print("\nConclusion: Since we have found an example of an infinite set (the positive x-axis) that is")
    print("transient, the statement that every transient set must be finite is false.")

solve_task()