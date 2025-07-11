import math

def calculate_average_projection_area(n):
    """
    Calculates the average projected area of an n-dimensional cross-polytope for odd n.
    """
    if not isinstance(n, int) or n < 3 or n % 2 == 0:
        print("Error: The dimension n must be an odd integer greater than or equal to 3.")
        return

    # k = (n-1)/2
    k = (n - 1) // 2

    # Numerator of the formula: n * sqrt(n)
    numerator = n * math.sqrt(n)
    
    # Denominator of the formula: 2 * (k!)^2
    try:
        k_factorial = math.factorial(k)
        denominator = 2 * (k_factorial**2)
    except ValueError:
        print(f"Error: k = {k} is too large for factorial calculation.")
        return
    
    # The final result
    average_area = numerator / denominator

    # Output the steps as requested
    print(f"For dimension n = {n}:")
    print(f"We use the formula: E[A] = (n * sqrt(n)) / (2 * (((n-1)/2)!)^2)")
    print(f"Here, k = (n-1)/2 = ({n}-1)/2 = {k}")
    print(f"The equation becomes: E[A] = ({n} * sqrt({n})) / (2 * ({k}!)^2)")
    print("-" * 20)
    print(f"Numerator = {n} * {math.sqrt(n):.6f} = {numerator:.6f}")
    print(f"Denominator = 2 * ({k_factorial})^2 = 2 * {k_factorial**2} = {denominator}")
    print(f"Average projected area = {numerator:.6f} / {denominator} = {average_area:.6f}")

if __name__ == '__main__':
    # Set the odd dimension n here. For example, n=3 for the octahedron.
    # The dimension must be an odd integer >= 3.
    n = 5
    calculate_average_projection_area(n)
