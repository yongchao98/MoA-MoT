import math

def solve():
    """
    This function formulates and prints the d-threshold for Hamiltonicity
    for the given range of d.
    """

    # The problem asks for the d-threshold for Hamiltonicity, p.
    # The value of d is given in terms of n and eta.
    # d = n/2 - eta

    # According to results in random graph theory, this threshold is governed by
    # the threshold for a random graph G(d, p) to have a Hamiltonian path.
    # This threshold p is given by the formula:
    # p = (ln(d) + ln(ln(d))) / d

    # We will format this equation as a string for the final output.
    # The instruction asks to output each number, which we interpret as
    # making sure all constants are clearly written in the final equation.
    
    # Define the components of the equation as strings
    var_p = "p"
    var_d = "d"
    var_n = "n"
    var_eta = "\\eta" # Using a LaTeX-style representation for eta
    number_2_str = "2"

    # Assemble the definition of d
    d_definition = f"{var_d} = {var_n}/{number_2_str} - {var_eta}"

    # Assemble the equation for p in terms of d
    p_equation_d = f"({var_p}) = (ln({var_d}) + ln(ln({var_d}))) / {var_d}"

    # Assemble the full equation for p in terms of n and eta
    d_expr = f"({var_n}/{number_2_str} - {var_eta})"
    numerator = f"ln({d_expr}) + ln(ln({d_expr}))"
    denominator = d_expr
    p_equation_full = f"{var_p} = ({numerator}) / {denominator}"
    
    print("The d-threshold for Hamiltonicity is given by the asymptotic formula below.")
    print("First, the variable d is defined as:")
    print(d_definition)
    print("\nThe threshold p is given in terms of d as:")
    print(p_equation_d)
    print("\nSubstituting d, the full expression for the threshold p is:")
    print(p_equation_full)

solve()