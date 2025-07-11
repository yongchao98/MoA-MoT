def generate_formula_string():
    """
    Constructs and prints the symbolic expression for the infinite product
    \prod_{n=3}^{\infty}\left(1-\frac{z^3}{n^3}\right).
    
    The code demonstrates how to build the string for the final formula,
    explicitly including all the numbers from the problem statement as requested.
    """
    
    # Define string representations for the variable and symbols
    z_str = "z"
    pi_str = "pi"
    I_str = "I"
    
    # Left hand side of the equation as given in the problem
    lhs_str = "prod_{n=3 to inf}(1 - z**3/n**3)"

    # Right hand side components
    # The terms for n=1 and n=2 that are divided out
    term_n1_str = f"(1 - {z_str}**3 / 1**3)"
    term_n2_str = f"(1 - {z_str}**3 / 2**3)"

    # The Gamma function product from the general identity
    gamma1_str = f"Gamma(1 - {z_str})"
    gamma2_str = f"Gamma(1 - {z_str}*exp(2*{pi_str}*{I_str}/3))"
    gamma3_str = f"Gamma(1 - {z_str}*exp(4*{pi_str}*{I_str}/3))"
    
    # Full right hand side of the equation
    rhs_str = f"1 / ({term_n1_str} * {term_n2_str} * {gamma1_str} * {gamma2_str} * {gamma3_str})"
    
    # Print the final equation
    print(f"{lhs_str} = {rhs_str}")

if __name__ == '__main__':
    generate_formula_string()