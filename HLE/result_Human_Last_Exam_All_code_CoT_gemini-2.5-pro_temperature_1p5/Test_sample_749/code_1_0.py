import math

def solve_branching_walk_limit():
    """
    Calculates the limiting return probability for the branching random walk.

    The probability of infinitely many particles visiting site 0 is non-zero only
    if the associated process is supercritical. We analyze the mean of this process
    in the limit as h->0.
    """
    
    # In the h->0 limit, the environment is all-blue.
    # Jump probabilities for a blue site:
    p_right = 4/5
    p_left = 1/5

    print("Step 1: Define jump probabilities in the h->0 (all-blue) limit.")
    print(f"Probability of jumping right (p): {p_right}")
    print(f"Probability of jumping left (q): {p_left}\n")

    # Let R_i be the probability that a walk starting at site i ever reaches site 0.
    # The return probability from site 0, p_return, is p * R_1 + q * R_{-1}.
    # The jump probabilities p and q from site 0 are p_right and p_left.
    p0 = p_right
    q0 = p_left

    # Step 2: Calculate R_1, the probability of reaching 0 from 1.
    # For a 1D random walk with p > q, the probability of ever hitting i-1 from i
    # is q/p. The drift is to the right, away from 0.
    R1 = p_left / p_right
    print("Step 2: Calculate the probability of reaching site 0 from site 1 (R_1).")
    print(f"The walk has a drift to the right, so R_1 = q / p.")
    print(f"R_1 = {p_left} / {p_right} = {R1}\n")

    # Step 3: Calculate R_{-1}, the probability of reaching 0 from -1.
    # The drift is to the right, towards 0. Hitting is certain.
    R_minus_1 = 1.0
    print("Step 3: Calculate the probability of reaching site 0 from site -1 (R_-1).")
    print("The drift is to the right, so the walk is certain to hit any site to its right.")
    print(f"R_-1 = {R_minus_1}\n")

    # Step 4: Calculate the return probability p_return for the limiting process.
    p_return = p0 * R1 + q0 * R_minus_1
    print("Step 4: Calculate the overall return probability to site 0, p_return.")
    print("p_return = p * R_1 + q * R_-1")
    # Using print to show the equation with numbers
    print(f"p_return = {p0} * {R1} + {q0} * {R_minus_1} = {p0 * R1 + q0 * R_minus_1}")

    print("\n--- Conclusion ---")
    print(f"The limiting mean number of returning offspring is {p_return:.2f}, which is less than 1.")
    print("This means the process of visiting site 0 is subcritical in the limit.")
    print("For a subcritical process, the probability of infinite survival (in this case, infinite visits) is 0.")

solve_branching_walk_limit()