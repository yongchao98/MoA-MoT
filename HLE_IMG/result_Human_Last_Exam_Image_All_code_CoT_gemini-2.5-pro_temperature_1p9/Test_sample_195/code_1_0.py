def generate_equation():
    """
    This function formulates and prints the equation for f(x) based on graphical analysis.
    """
    # Numerator derived from x-intercepts (-b, a, b)
    # The intercept at 'a' is a turning point, so it's a double root (x-a)^2.
    # The intercepts at '-b' and 'b' are simple crossings, so (x+b)(x-b) = x^2 - b^2.
    numerator_str = "(x - a)^2 * (x^2 - b^2)"

    # Denominator derived from vertical asymptotes (c, d)
    # At x=d, the function goes to +infinity on both sides, suggesting an even power (x-d)^2.
    # At x=c, the function goes to +infinity on one side and -infinity on the other, suggesting an odd power (x-c).
    denominator_str = "(x - d)^2 * (x - c)"

    # Determine the length for the fraction line for neat printing
    line_length = max(len(numerator_str), len(denominator_str))

    # Create the fraction line string
    fraction_line = "-" * (line_length + 2)

    # Center the numerator and denominator strings for better alignment
    centered_numerator = numerator_str.center(line_length)
    centered_denominator = denominator_str.center(line_length)

    # Print the final equation in a structured format
    print("Based on the graphical features, a possible equation for the function f(x) is:")
    print()
    print(f"         {centered_numerator}")
    print(f"f(x) = {fraction_line}")
    print(f"         {centered_denominator}")

generate_equation()