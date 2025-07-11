def solve_random_walk_limit():
    """
    This function calculates the limit of the asymptotic speed v(c) as c -> infinity.
    
    The final speed is given by the formula:
    v = pi_0 * v_0 + pi_1 * v_1
    where:
    - v_0 is the speed on the lower rail.
    - v_1 is the speed on the upper rail.
    - pi_0 and pi_1 are the stationary probabilities of being on the lower/upper rail.
    """

    # As c -> infinity, the speed on the lower rail approaches 1.
    v_0 = 1
    
    # As c -> infinity, the speed on the upper rail is 0 due to trapping.
    v_1 = 0
    
    # As c -> infinity, the stationary probability of being on the lower rail approaches 1.
    pi_0 = 1
    
    # As c -> infinity, the stationary probability of being on the upper rail approaches 0.
    pi_1 = 0
    
    # The final limiting speed is the weighted average.
    limit_v = pi_0 * v_0 + pi_1 * v_1
    
    print("Finding the limit of the asymptotic speed v(c) as c -> infinity:")
    print("-" * 60)
    print("The final speed is calculated as a weighted average of the speeds on each rail.")
    print("v = pi_0 * v_0 + pi_1 * v_1\n")
    
    print(f"The limiting speed on the lower rail (v_0) is: {v_0}")
    print(f"The limiting speed on the upper rail (v_1) is: {v_1}")
    print(f"The limiting probability of being on the lower rail (pi_0) is: {pi_0}")
    print(f"The limiting probability of being on the upper rail (pi_1) is: {pi_1}\n")

    print("Substituting these values into the equation for the limit:")
    print(f"lim v(c) = ({pi_0} * {v_0}) + ({pi_1} * {v_1})")
    print(f"lim v(c) = {pi_0 * v_0} + {pi_1 * v_1}")
    print(f"lim v(c) = {limit_v}")
    
solve_random_walk_limit()