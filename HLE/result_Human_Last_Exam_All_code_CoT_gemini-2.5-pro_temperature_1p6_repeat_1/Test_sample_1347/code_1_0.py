def display_r0f_expression():
    """
    This function formats and prints the derived expression for R0f.
    Since no numerical values are given, it uses symbolic variable names.
    """
    
    # Define symbolic representations of the variables as strings
    b = "b"
    pg = "pg"
    c = "c"
    pt = "pt"
    tau = "𝜏"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    mu_g = "𝜇g"
    
    # Construct the numerator and denominator of the expression
    # The order of terms in the numerator can vary, but the result is the same.
    # We follow the logic of the derivation: (b*pg) * (c*pt) * tau
    numerator = f"{b} * {pg} * {c} * {pt} * {tau}"
    
    # Denominator terms correspond to the rates of leaving each state
    denominator_part1 = f"({gamma_t} + {mu_t})"
    denominator_part2 = f"({tau} + {mu_g})"
    denominator_part3 = f"{mu_g}"
    denominator = f"{denominator_part1} * {denominator_part2} * {denominator_part3}"
    
    # Print the final expression in a clear format
    print("The derived expression for R0f is:")
    print()
    print(f"      {numerator}")
    print(f"R0f = —————————————————————————————————————————————————————")
    print(f"      {denominator}")

if __name__ == "__main__":
    display_r0f_expression()