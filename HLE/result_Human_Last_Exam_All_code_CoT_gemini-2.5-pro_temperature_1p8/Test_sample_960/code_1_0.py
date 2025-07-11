import sympy

def solve_probability():
    """
    Calculates the probability that the number of good products will equal
    the number of defective products.
    """
    v = sympy.Symbol('v')

    # The limiting proportion of good products follows a Beta(2, 1) distribution.
    # The PDF of a Beta(2, 1) distribution is f(v) = 2v for v in [0, 1].
    pdf = 2 * v

    # The conditional probability of the numbers becoming equal (P(T < ∞ | V=v)) is:
    # 1 if v <= 1/2
    # (1-v)/v if v > 1/2
    prob_cond_part1 = 1
    prob_cond_part2 = (1 - v) / v

    # We calculate the total probability by integrating P(T < ∞ | V=v) * f(v) dv
    # over the range of v.
    
    # Part 1: Integral from v = 0 to 1/2
    integrand1 = prob_cond_part1 * pdf
    integral1 = sympy.integrate(integrand1, (v, 0, 1/2))

    # Part 2: Integral from v = 1/2 to 1
    integrand2 = prob_cond_part2 * pdf
    integral2 = sympy.integrate(integrand2, (v, 1/2, 1))

    # The total probability is the sum of the two integrals.
    total_probability = integral1 + integral2

    print("This problem is equivalent to a Polya's Urn problem.")
    print("We start with W_0 = 2 good and B_0 = 1 defective products.")
    print("The goal is to find the upper bound for the probability of reaching a state where W_t = B_t.")
    print("\nWe calculate the exact probability, which is the least upper bound.")
    print("The calculation involves integrating a conditional probability over the Beta(2, 1) distribution of the limiting proportion 'v'.")
    print("\nThe probability P is the sum of two integrals:")
    print("P = Integral from 0 to 1/2 of [1 * (2v)] dv + Integral from 1/2 to 1 of [(1-v)/v * (2v)] dv")
    
    print("\nStep 1: Calculate the first integral.")
    print(f"Integral(2v, (v, 0, 1/2)) = [v^2] from 0 to 1/2 = (1/2)^2 - 0^2")
    
    print("\nStep 2: Calculate the second integral.")
    print(f"Integral(2(1-v), (v, 1/2, 1)) = [2v - v^2] from 1/2 to 1 = (2-1) - (1-1/4)")
    
    print("\nStep 3: Sum the results to get the final probability.")
    # We output each number in the final equation as requested.
    print(f"P = {integral1} + {integral2} = {total_probability}")
    
    print(f"\nThe exact probability is {total_probability}. This value is the least upper bound.")
    print(f"As a decimal, this is {float(total_probability)}.")

solve_probability()
<<<0.5>>>