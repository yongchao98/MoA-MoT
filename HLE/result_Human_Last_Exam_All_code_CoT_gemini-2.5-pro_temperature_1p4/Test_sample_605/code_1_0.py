import math
from functools import reduce

def solve_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant, interpreted as the Milnor number,
    for the given Calabi-Yau link.
    """
    # The weights and polynomial are given in the problem description.
    # We will use the weights to determine the degree 'd' based on the Calabi-Yau condition.
    weights = [22, 29, 49, 50, 75]
    d = sum(weights)

    # The Milnor number formula for a quasi-homogeneous polynomial singularity is:
    # μ = product over i of (d/w_i - 1)
    
    numerators = [d - w for w in weights]
    denominators = weights

    # We calculate the final result by multiplying the fractions.
    # Using integer arithmetic throughout avoids floating point errors.
    numerator_product = reduce(lambda x, y: x*y, numerators)
    denominator_product = reduce(lambda x, y: x*y, denominators)
    result = numerator_product // denominator_product

    # The final output will be formatted to show the calculation steps.
    print("The Crawley-Nordström invariant for the given Calabi-Yau link is interpreted as the Milnor number (μ).")
    print("The calculation is based on the weights and the quasi-homogeneous degree of the polynomial.\n")
    print(f"Weights: w = {tuple(weights)}")
    print(f"Degree: d = {' + '.join(map(str, weights))} = {d}\n")
    print("The Milnor number is calculated using the formula:")
    print("μ = (d/w₁ - 1) * (d/w₂ - 1) * (d/w₃ - 1) * (d/w₄ - 1) * (d/w₅ - 1)\n")
    
    print("Substituting the values into the formula:")
    
    # Building the string for the equation with fractions
    frac_terms = [f"(({d} - {w})/{w})" for w in weights]
    equation_str = " * ".join(frac_terms)
    print(f"μ = {equation_str}")
    
    final_equation_terms = [f"({num}/{den})" for num, den in zip(numerators, denominators)]
    final_equation_str = " * ".join(final_equation_terms)
    
    print(f"μ = {final_equation_str}")
    
    # The final equation as requested, showing each number.
    print(f"μ = {numerators[0]}/{denominators[0]} * {numerators[1]}/{denominators[1]} * {numerators[2]}/{denominators[2]} * {numerators[3]}/{denominators[3]} * {numerators[4]}/{denominators[4]} = {result}")

solve_crawley_nordstrom_invariant()
<<<1568>>>