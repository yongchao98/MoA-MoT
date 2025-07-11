import sympy

def solve_branching_walk():
    """
    Solves the branching random walk problem by calculating the criticality parameter.
    """
    # Define h as a symbol for our calculations
    h = sympy.Symbol('h')

    # Step 1: Define the mean number of offspring, m(h)
    # A particle becomes 1 with probability (1-h) and 2 with probability h.
    mean_offspring = 1 * (1 - h) + 2 * h
    
    # Step 2: Define the averaged jump probabilities, p_L(h) and p_R(h)
    # p_L = P(Left|Red)P(Red) + P(Left|Blue)P(Blue)
    prob_left = sympy.Rational(4, 5) * h + sympy.Rational(1, 5) * (1 - h)
    
    # p_R = P(Right|Red)P(Red) + P(Right|Blue)P(Blue)
    prob_right = sympy.Rational(1, 5) * h + sympy.Rational(4, 5) * (1 - h)

    # Step 3: Define the probability decay factor, gamma(h)
    gamma = 2 * sympy.sqrt(prob_left * prob_right)

    # Step 4: Define the criticality parameter, K(h) = m(h) * gamma(h)
    criticality_param = mean_offspring * gamma

    # Step 5: Evaluate K(h) in the limit h -> 0 by substituting h=0
    val_at_zero = criticality_param.subs(h, 0)

    # Step 6: Print the reasoning based on the result
    print("Let P(E) be the probability that site 0 is visited by infinitely many particles.")
    print("We want to calculate lim_{h->0} P(E).")
    print("\nFor a branching random walk, a site is visited infinitely often only if a criticality parameter K(h) >= 1.")
    print("If K(h) < 1, the site is visited only a finite number of times, and thus P(E) = 0.")
    print("\nThis parameter K(h) is the product of two terms:")
    print("1. The mean number of offspring: m(h) = (1-h)*1 + h*2")
    print("2. The single-walker decay factor: gamma(h) = 2 * sqrt(p_L(h) * p_R(h))")
    
    print("\nLet's evaluate K(h) for h -> 0 by setting h=0.")

    # Get the numeric values for the final equation
    m_0 = mean_offspring.subs(h,0)
    pl_0 = prob_left.subs(h,0)
    pr_0 = prob_right.subs(h,0)
    gamma_0 = gamma.subs(h,0)

    print("\nFinal calculation:")
    print(f"At h=0, the mean offspring m(0) = {m_0}.")
    print(f"The left-jump probability p_L(0) = {pl_0}.")
    print(f"The right-jump probability p_R(0) = {pr_0}.")
    print(f"The decay factor gamma(0) = 2 * sqrt({pl_0} * {pr_0}) = {gamma_0}.")
    
    print(f"\nThe criticality parameter at h=0 is:")
    print(f"K(0) = m(0) * gamma(0)")
    print(f"K(0) = {m_0} * {gamma_0} = {val_at_zero}")

    print(f"\nSince K(0) = {val_at_zero} < 1, K(h) remains less than 1 for all sufficiently small h > 0.")
    print("This implies that for small h, P(E) = 0.")
    print("Therefore, the limit of P(E) as h approaches 0 is 0.")

solve_branching_walk()
<<<0>>>