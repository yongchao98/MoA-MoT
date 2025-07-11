from fractions import Fraction

def solve_speed():
    """
    Calculates the asymptotic speed v(c) as c -> infinity.
    """
    # Probabilities of edges being present
    # p_R: Probability of an upper horizontal edge existing
    p_R = Fraction(2, 3)
    # p_V: Probability of a vertical edge existing
    p_V = Fraction(1, 2)

    # Step 1: Calculate expected time t_0 to cross one column from the lower rail.
    # From (n,0), the edge to (n+1,0) always exists. The weight e^c for this
    # move is infinitely larger than for any other move. So the walker moves
    # (n,0)->(n+1,0) in one step.
    t_0 = Fraction(1)

    # Step 2: Calculate expected time t_1 to cross one column from the upper rail.
    # Let t_1 be the expected time starting from (n,1) to reach column n+1.
    # At (n,1), the walker's choice depends on the available edges:
    # Case 1: Right edge exists (prob p_R = 2/3). Walker takes it. Time = 1.
    # Case 2: Right edge absent, Down edge exists (prob (1-p_R)*p_V = 1/3 * 1/2 = 1/6).
    #         Walker goes (n,1)->(n,0), then (n,0)->(n+1,0). Time = 1 + t_0 = 2.
    # Case 3: Right and Down edges absent (prob (1-p_R)*(1-p_V) = 1/3 * 1/2 = 1/6).
    #         Walker must go left to (n-1,1). Time = 1 + E[time from (n-1,1) to n+1].
    # The time from (n-1,1) to n+1 is t_1 + E[time from column n].
    # To find E[time from column n], we need the probability of landing at (n,1) vs (n,0)
    # when starting from (n-1,1).
    # Let pi_1 be the probability of landing on the upper rail.
    # pi_1 = p_R (go right) + (1-p_R)*(1-p_V)*pi_1 (go left, then succeed from there).
    # pi_1 = 2/3 + 1/6 * pi_1 => (5/6)pi_1 = 2/3 => pi_1 = 4/5.
    pi_1_from_left = Fraction(4, 5)
    pi_0_from_left = 1 - pi_1_from_left
    #
    # The recurrence for t_1 is:
    # t_1 = p_R*1 + (1-p_R)*p_V*(1+t_0) + (1-p_R)*(1-p_V)*(1 + t_1 + pi_1_from_left*t_1 + pi_0_from_left*t_0)
    # t_1 = 2/3*1 + 1/6*(2) + 1/6*(1 + t_1 + 4/5*t_1 + 1/5*1)
    # t_1 = 1 + 1/6 * (6/5 + 9/5 * t_1)
    # t_1 = 1 + 1/5 + 3/10 * t_1
    # t_1 * (7/10) = 6/5
    t_1 = Fraction(6, 5) / Fraction(7, 10)

    # Step 3: Define the Markov Chain for rail transitions and find stationary distribution.
    # p_ij = prob of entering col n+1 on rail j, given you entered col n on rail i.
    # From rail 0, you always stay on rail 0.
    p_00, p_01 = Fraction(1), Fraction(0)
    # From rail 1, we need to find the probability of landing on rail 1 next (p_11).
    # This calculation is similar to finding pi_1_from_left.
    # p_11 = p_R*1 + (1-p_R)*(1-p_V)*(pi_1_from_left*p_11 + pi_0_from_left*p_01)
    # p_11 = 2/3 + 1/6 * (4/5 * p_11 + 1/5 * 0)
    # p_11 = 2/3 + 2/15 * p_11 => (13/15)p_11 = 2/3
    p_11 = Fraction(2, 3) / Fraction(13, 15)
    p_10 = 1 - p_11
    
    # The transition matrix is [[1, 0], [p_10, p_11]]. State 0 is absorbing.
    # The stationary distribution (pi_0, pi_1) must satisfy pi*P = pi.
    # pi_0 = pi_0*p_00 + pi_1*p_10 => pi_0 = pi_0*1 + pi_1*p_10 => pi_1*p_10 = 0
    # Since p_10 = 3/13 != 0, this implies pi_1 = 0.
    pi_1 = Fraction(0)
    # As pi_0 + pi_1 = 1, we have pi_0 = 1.
    pi_0 = Fraction(1)

    # Step 4: Calculate the final speed v = 1 / (pi_0*t_0 + pi_1*t_1)
    avg_time_per_col = pi_0 * t_0 + pi_1 * t_1
    speed = 1 / avg_time_per_col
    
    print("This script calculates the asymptotic speed of the random walk.")
    print("The final speed v is given by the formula v = 1 / (pi_0*t_0 + pi_1*t_1)")
    print("\nCalculated values for each number in the final equation:")
    print(f"t_0 (Expected time from lower rail): {t_0}")
    print(f"t_1 (Expected time from upper rail): {t_1}")
    print(f"pi_0 (Stationary probability for lower rail): {pi_0}")
    print(f"pi_1 (Stationary probability for upper rail): {pi_1}")
    print("\nFinal calculation:")
    print(f"v = 1 / ({pi_0}*{t_0} + {pi_1}*{t_1})")
    print(f"v = 1 / {avg_time_per_col}")
    print(f"v = {speed}")

solve_speed()
<<<1>>>