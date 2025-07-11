import math

def solve():
    """
    Solves the branching random walk problem.

    The problem asks for the limit as h->0 of the probability that site 0 is
    visited by infinitely many particles.

    Let P(h) be this probability for a given h.

    1.  At each site i, the local drift is E[X_{n+1}-X_n | X_n=i].
        - If site i is blue (prob 1-h), drift = (1)*(4/5) + (-1)*(1/5) = 3/5.
        - If site i is red (prob h), drift = (1)*(1/5) + (-1)*(4/5) = -3/5.

    2.  The average, or "annealed", drift over the random environment is:
        E[drift] = (1-h) * (3/5) + h * (-3/5) = (3/5) * (1 - 2*h).

    3.  For a one-dimensional branching random walk in a random environment, the
        asymptotic velocity of the particle cloud has the same sign as the
        annealed drift. Let v_min(h) be the velocity of the leftmost particle.
        If E[drift] > 0, then v_min(h) > 0.

    4.  The condition E[drift] > 0 translates to (3/5)*(1-2*h) > 0, which
        is true for h < 1/2.

    5.  The limit h->0 considers h in the region (0, 1/2). For any such h,
        the velocity of the leftmost particle is positive (v_min(h) > 0). This
        means the entire cloud of particles almost surely drifts to positive
        infinity.

    6.  If the entire population drifts to +infinity, any fixed site (like site 0)
        will be visited only a finite number of times. The probability of
        infinite visits is therefore 0.

    7.  So, for any h in (0, 1/2), P(h) = 0.

    8.  The limit is lim_{h->0} P(h) = lim_{h->0} 0 = 0.
    """

    # The final probability is P.
    # The governing logic leads to the conclusion that P=0.
    # We will print the numbers in this final equation.
    final_probability = 0
    
    print("The final probability is P, determined by the asymptotic behavior of the system.")
    print("The analysis shows that for any h in (0, 1/2), the particle population drifts to +infinity.")
    print("This means the probability of infinitely many visits to site 0 is 0.")
    print("The limiting probability as h->0 is therefore 0.")
    print("Final equation: P = 0")
    print("The number in the final equation is:")
    print(final_probability)

solve()