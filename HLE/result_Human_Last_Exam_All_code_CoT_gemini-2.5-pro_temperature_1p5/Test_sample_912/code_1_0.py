def print_work_formula():
    """
    Prints the symbolic formula for the work done by the current source in the given cycle.
    """
    
    # Define symbolic components of the formula
    numerator_factor1 = "(μ - μ_0)"
    numerator_factor2 = "N²"
    numerator_factor3 = "w"
    numerator_factor4 = "(x₂ - x₁)"
    numerator_factor5 = "(I₂² - I₁²)"
    
    denominator = "2g"
    
    # Assemble the final formula string
    work_formula = (
        f"W = - ( {numerator_factor1} / {denominator} ) * "
        f"{numerator_factor2} * {numerator_factor3} * {numerator_factor4} * {numerator_factor5}"
    )

    # Reformat for a cleaner fractional representation
    final_formula = (
        f"W = - (({numerator_factor1}) * {numerator_factor2} * {numerator_factor3} * {numerator_factor4} * {numerator_factor5}) / ({denominator})"
    )
    
    print("The work done by the current source for each cycle is given by the formula:")
    # A slightly more readable version:
    print(f"W = - \u0332{numerator_factor1}\u0332 \u00B7 N\u00B2 \u00B7 w \u00B7 (x\u2082 - x\u2081) \u00B7 (I\u2082\u00B2 - I\u2081\u00B2)")
    print(f"      {denominator.center(len(numerator_factor1) + 2)}")

    # And now presenting the final answer in the single-line format for clarity, matching choice D
    final_expression = f"W = - ((\u03BC - \u03BC\u2080) / (2*g)) * N**2 * w * (x\u2082 - x\u2081) * (I\u2082**2 - I\u2081**2)"
    
    print("\nOr, written in a single line:")
    # This simplified print statement matches option D directly.
    # The term (\mu - \mu_0)/(2g) is a single fraction.
    print("W = - ((\u03BC - \u03BC\u2080) * N\u00B2 * w * (x\u2082 - x\u2081) * (I\u2082\u00B2 - I\u2081\u00B2)) / (2g)")


if __name__ == "__main__":
    print_work_formula()