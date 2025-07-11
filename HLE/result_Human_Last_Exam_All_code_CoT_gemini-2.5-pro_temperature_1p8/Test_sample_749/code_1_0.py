import math

def solve_branching_walk():
    """
    Calculates the limit of the probability that site 0 is visited by
    infinitely many particles in a branching random walk.
    """
    # In the limit h -> 0, the environment is all blue.
    p_left_blue = 1/5
    p_right_blue = 4/5

    # 1. Calculate the drift of the random walk.
    v = p_right_blue - p_left_blue

    # 2. Calculate the decay factor for the probability of visiting site 0.
    # The probability of visiting 0 from site x is (p_left/p_right)^x = (1/k)^x.
    k = p_right_blue / p_left_blue

    # 3. Check the convergence condition for the expected number of visitors.
    # The expected number is finite if 1 < k^v (for h -> 0).
    k_to_the_v = k**v

    print("Step 1: In the h->0 limit, we analyze the system in an all-blue environment.")
    print(f"The probability of moving left is p = {p_left_blue}")
    print(f"The probability of moving right is q = {p_right_blue}")
    print("-" * 20)
    
    print("Step 2: Calculate the drift v of the particle cloud.")
    print(f"v = q - p = {p_right_blue} - {p_left_blue} = {v}")
    print("-" * 20)

    print("Step 3: Determine the probability of visiting site 0.")
    print("The probability for a particle at site x to visit site 0 is (p/q)^x = (1/k)^x.")
    print(f"k = q / p = {p_right_blue} / {p_left_blue} = {k}")
    print("-" * 20)
    
    print("Step 4: Check the condition for the expected number of visitors to be finite.")
    print("The expectation is finite if 1 + h < k^v. As h -> 0, this is 1 < k^v.")
    print("The final equation to check is:")
    print(f"1 < {k}^({v})")
    print("-" * 20)
    
    print("Step 5: Evaluating the terms in the equation.")
    print(f"{k}^({v}) = {k_to_the_v:.4f}")
    
    is_condition_met = 1 < k_to_the_v
    print(f"Is 1 < {k_to_the_v:.4f}? {is_condition_met}")
    print("-" * 20)
    
    print("Conclusion: Since the condition is true, the expected number of particles visiting site 0 is finite.")
    print("By the first Borel-Cantelli lemma, if the expected number of events is finite, the probability of infinitely many events occurring is 0.")
    
    final_answer = 0
    print(f"\nFinal calculated limit: {final_answer}")

solve_branching_walk()