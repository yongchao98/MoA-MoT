import math

def solve_dna_probability():
    """
    Calculates the limiting probability for the quantum DNA polymerase problem.
    """
    # The system has M distinct bases.
    M = 8

    # The problem asks for the limiting value of P(n), the probability that for a random
    # DNA sequence of length n, the condition (count of base x) mod M != x for all x.
    
    # As derived in the thinking steps, the exact probability P(n) oscillates slightly
    # based on the value of n mod M. The problem asks for a single limiting value,
    # which is best interpreted as the average probability over these oscillations.
    
    # This average limiting probability simplifies to the elegant closed-form expression ((M-1)/M)^M.

    # Numerator of the base of the expression
    num_base = M - 1
    # Denominator of the base of the expression
    den_base = M
    # The exponent in the expression
    power = M

    # Calculate the numerator and denominator of the final probability
    final_numerator = num_base ** power
    final_denominator = den_base ** power

    # Calculate the final probability as a decimal
    probability = final_numerator / final_denominator

    # Output the steps of the calculation as requested
    print(f"The number of bases is M = {M}.")
    print(f"The limiting average probability P is given by the expression ( (M-1) / M )^M.")
    print(f"Substituting M = {M}, we get P = ({num_base}/{den_base})^{power}.")
    print(f"The final equation is P = {num_base}^{power} / {den_base}^{power}.")
    print(f"P = {final_numerator} / {final_denominator}")
    print(f"The numerical value is P \u2248 {probability}")

solve_dna_probability()