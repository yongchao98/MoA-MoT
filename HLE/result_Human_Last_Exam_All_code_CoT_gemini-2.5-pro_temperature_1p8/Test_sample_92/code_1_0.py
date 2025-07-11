import sympy

def solve():
    """
    Calculates the probability that the marble escapes.
    
    The problem can be solved by considering the ratio of probabilities of jumping
    to the escape bin versus the melt bin in a single step from any position n < 2024.
    
    Let n be the current position of the marble (n < 2024).
    The escape bin is at N_E = 2025.
    The melting bin is at N_M = 2024.
    
    The probability of jumping by i is (1/3)^|i|.
    
    Prob(jump to N_E from n) = (1/3)^(N_E - n)
    Prob(jump to N_M from n) = (1/3)^(N_M - n)
    
    The conditional probability of escaping, given that the process terminates
    in the next step, is:
    P(escape | absorbed) = Prob(jump to N_E) / (Prob(jump to N_E) + Prob(jump to N_M))
                        = (1/3)^(2025-n) / ((1/3)^(2025-n) + (1/3)^(2024-n))
                        = (1/3) / (1/3 + 1)
                        = (1/3) / (4/3)
                        = 1/4
    
    This ratio is independent of the starting position n. This implies that
    at any point in its random walk, the marble's odds of escaping vs. melting,
    should the process end, are the same. Thus, the overall probability is this ratio.
    """
    
    # Let's represent the probabilities as fractions for exact arithmetic.
    # prob_jump_right = 1/3 (this is the extra factor for 2025 vs 2024)
    # prob_jump_left = 1
    prob_to_escape_factor = sympy.Rational(1, 3)
    prob_to_melt_factor = sympy.Rational(1, 1)

    # The probability of escape is the ratio of the escape factor to the sum of both factors.
    total_factor = prob_to_escape_factor + prob_to_melt_factor
    probability = prob_to_escape_factor / total_factor

    print(f"The numerator in the final fraction is the relative probability factor of jumping to the escape portal compared to the melting torch.")
    print(f"Let the base probability of reaching the torch at bin 2024 be P. The probability of reaching the escape portal at bin 2025 is P * (1/3).")
    print(f"The probability of escaping is P*(1/3) / (P + P*(1/3)) which simplifies to (1/3) / (1 + 1/3).")

    num = prob_to_escape_factor.p
    den_part1 = prob_to_melt_factor.p
    den_part2 = prob_to_escape_factor.p

    print(f"The final simplified equation is: {num} / ({den_part1} + {num}/3 * 3) = {num} / ({den_part1} + {den_part2}/3*3) which is {prob_to_escape_factor} / {total_factor}")
    print(f"P(escape) = {prob_to_escape_factor} / ({prob_to_melt_factor} + {prob_to_escape_factor}) = {probability}")

solve()
<<<0.25>>>