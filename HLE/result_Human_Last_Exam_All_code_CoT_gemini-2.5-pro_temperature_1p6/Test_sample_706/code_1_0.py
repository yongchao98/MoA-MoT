import math

def solve():
    """
    This function explains and calculates the limit of the asymptotic speed v(c) as c approaches infinity.
    
    Let v_0 be the speed on the lower level and v_1 be the speed on the upper level.
    As c -> infinity, the walk is dominated by moves that maximize horizontal displacement.

    1. Speed on the lower level (v_0):
    A walker at (n, 0) can move to (n+1, 0) [weight e^c], (n-1, 0) [weight e^(-c)], or (n, 1) [weight 1].
    The probability of moving right P_right = e^c / (e^c + e^(-c) + 1) -> 1 as c -> infinity.
    Each step advances the position by 1. So, the speed v_0 = 1.

    2. Speed on the upper level (v_1):
    An excursion on the upper level starts and ends on the lower level. We can calculate the expected time (T) and displacement (X) for such an excursion before returning to the lower level.
    Using recurrence relations based on the probabilities of edges existing (p_H = 2/3 for horizontal, p_V = 1/2 for vertical):
    For c -> infinity, the greedy path is taken.
    T = (2/3)*(1+T) + (1/3)*[ (1/2)*(1) + (1/2)*(1+T) ]
    T = 2/3 + 2T/3 + 1/6 + 1/6 + T/6
    T = 1 + 5T/6 => T/6 = 1 => T = 6.
    
    X = (2/3)*(1+X) + (1/3)*[ (1/2)*(0) + (1/2)*(-1+X) ]
    X = 2/3 + 2X/3 - 1/6 + X/6
    X = 1/2 + 5X/6 => X/6 = 1/2 => X = 3.
    
    The speed on the upper level is v_1 = X/T = 3/6 = 1/2.

    3. Overall Speed v(c):
    The overall speed is the weighted average of the speeds on both levels, weighted by the stationary probability of being on each level (pi_0, pi_1).
    The transition probability from lower to upper is proportional to e^(-c), while the transition from upper to lower is a constant.
    Therefore, the stationary probability of being on the upper level, pi_1, approaches 0 as c -> infinity.
    pi_0 approaches 1.

    v(c) = pi_0 * v_0 + pi_1 * v_1
    lim_{c->inf} v(c) = 1 * v_0 + 0 * v_1 = 1 * 1 + 0 * 0.5 = 1.
    """
    
    # In the limit c -> infinity, the problem simplifies and can be solved analytically.
    
    # Speed on the lower level
    v_0 = 1
    
    # Expected time and displacement for an excursion on the upper level
    # From solving the recurrence T = 1 + (5/6)T
    T_upper_excursion = 6
    # From solving the recurrence X = 1/2 + (5/6)X
    X_upper_excursion = 3
    
    # Speed on the upper level
    v_1 = X_upper_excursion / T_upper_excursion
    
    # The stationary probability of being on the upper level (pi_1) approaches 0
    # as the probability of transitioning from lower to upper is proportional to exp(-c)
    pi_1_limit = 0
    
    # The stationary probability of being on the lower level (pi_0) approaches 1
    pi_0_limit = 1
    
    # The final asymptotic speed is the weighted average
    final_speed = pi_0_limit * v_0 + pi_1_limit * v_1
    
    print("Let T be the expected time and X be the expected horizontal displacement for an excursion on the upper level.")
    print(f"Solving the recurrence relations for the c -> infinity limit, we get:")
    print(f"E[T] = {T_upper_excursion}")
    print(f"E[X] = {X_upper_excursion}")
    print(f"The effective speed during an upper-level excursion is v1 = X/T = {X_upper_excursion}/{T_upper_excursion} = {v_1}")
    print("The speed on the lower level, v0, is 1.")
    print("The stationary probability of being on the upper level, pi_1, approaches 0 as c -> infinity.")
    print("The final asymptotic speed v = pi_0*v0 + pi_1*v1.")
    print(f"In the limit, this becomes 1 * {v_0} + 0 * {v_1} = {final_speed}")

solve()