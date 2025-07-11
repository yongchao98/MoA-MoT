import math

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau link.
    """
    # The weights of the ambient space P(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]
    n = len(weights)

    # The degree 'd' of the quasihomogeneous polynomial.
    # As determined by analyzing the monomials, d=225.
    d = 225

    # Check if the Calabi-Yau condition (d = sum of weights) is met
    if d != sum(weights):
        print(f"Warning: The Calabi-Yau condition is not met. d={d}, sum(w_i)={sum(weights)}")
    
    # Calculate the charges q_i = w_i / d
    charges = [w / d for w in weights]

    # Calculate the individual terms (1 - 2*q_i) for the invariant sum
    cn_terms = [1 - 2 * q for q in charges]

    # Sum the terms to get the final invariant
    cn_invariant = sum(cn_terms)

    # Build and print the equation as requested
    equation_parts = []
    for w in weights:
        equation_parts.append(f"(1 - 2 * {w}/{d})")
    
    equation_string = " + ".join(equation_parts)

    print("The Crawley-Nordström invariant is calculated using the formula:")
    print("c_N = sum_{i=1 to 5} (1 - 2 * w_i/d)\n")
    print("Substituting the given values:")
    # Using math.isclose() for a robust floating-point comparison to check if the result is an integer
    if math.isclose(cn_invariant, round(cn_invariant)):
        print(f"c_N = {equation_string} = {int(round(cn_invariant))}")
    else:
        print(f"c_N = {equation_string} = {cn_invariant}")

calculate_crawley_nordstrom_invariant()