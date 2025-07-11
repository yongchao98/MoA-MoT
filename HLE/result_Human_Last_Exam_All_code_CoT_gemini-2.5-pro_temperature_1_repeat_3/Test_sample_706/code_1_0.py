import math

def solve_speed_limit():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.
    """
    # --- Step 1: Define probabilities of edge existence ---
    p_H = 2/3  # Probability that an upper horizontal edge exists
    p_V = 1/2  # Probability that a vertical edge exists

    # --- Step 2: Calculate expected time T0 to cross one unit starting at level 0 ---
    # At (n,0), the edge to (n+1,0) always exists. As c->inf, the walker
    # takes this edge.
    T0 = 1.0

    # --- Step 3: Calculate expected time T1 to cross one unit starting at level 1 ---
    # We solve the equation for T1:
    # T1 = p_H * 1 (move right)
    #      + (1-p_H) * p_V * (1 + T0) (blocked right, move down then right)
    #      + (1-p_H) * (1-p_V) * (1 + 2*T1) (trapped, move left, then cross two units)
    # T1 = (2/3)*1 + (1/3)*(1/2)*(1+T0) + (1/3)*(1/2)*(1+2*T1)
    # T1 = 2/3 + 1/6 * (1+1) + 1/6 * (1+2*T1)
    # T1 = 2/3 + 1/3 + 1/6 + T1/3
    # T1 = 1 + 1/6 + T1/3
    # (2/3)*T1 = 7/6  => T1 = (7/6)*(3/2) = 7/4
    T1 = 7.0 / 4.0

    # --- Step 4: Determine the stationary distribution of arrival levels ---
    # Let pi_prime = (pi0, pi1) be the stationary probability of arriving at a
    # new column on level 0 or 1. We find this from the transition matrix M.
    # M[i][j] is the probability of arriving at level j starting from level i.
    # From level 0, the walker always moves to (n+1,0), so it arrives at level 0.
    M00, M01 = 1.0, 0.0
    
    # From level 1, the arrival level depends on the path taken.
    # This leads to a quadratic equation for M11: M11^2 - 6*M11 + 4 = 0
    # The solutions are 3 +/- sqrt(5). Only 3 - sqrt(5) is a valid probability.
    M11 = 3.0 - math.sqrt(5)
    M10 = 1.0 - M11
    
    # The stationary distribution (pi0, pi1) satisfies (pi0, pi1) * M = (pi0, pi1).
    # pi0 * M00 + pi1 * M10 = pi0  => pi0 + pi1*M10 = pi0 => pi1*M10 = 0.
    # Since M10 is not 0, pi1 must be 0.
    pi1_prime = 0.0
    pi0_prime = 1.0 - pi1_prime
    
    # --- Step 5: Calculate the average time per unit and the final speed ---
    # The average time is T_avg = pi0_prime * T0 + pi1_prime * T1
    T_avg = pi0_prime * T0 + pi1_prime * T1
    
    # The speed is the reciprocal of the average time.
    v_limit = 1.0 / T_avg
    
    # --- Step 6: Print the results ---
    print("This script calculates the asymptotic speed v of a biased random walk on a random ladder graph in the limit c -> infinity.")
    print("\n1. Expected time to cross one horizontal unit:")
    print(f"   - Starting on lower level (T0): {T0}")
    print(f"   - Starting on upper level (T1): {T1}")

    print("\n2. Stationary probabilities of arriving at a column on a given level:")
    print(f"   - Probability of arriving on lower level (pi'_0): {pi0_prime}")
    print(f"   - Probability of arriving on upper level (pi'_1): {pi1_prime}")

    print("\n3. Final speed calculation:")
    print(f"   v = 1 / (pi'_0 * T0 + pi'_1 * T1)")
    print(f"   v = 1 / ({pi0_prime} * {T0} + {pi1_prime} * {T1})")
    print(f"   v = 1 / {T_avg}")
    print(f"   v = {v_limit}")

solve_speed_limit()
<<<1.0>>>