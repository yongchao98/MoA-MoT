import math

def calculate_average_projection_area():
    """
    Calculates the average area of a projection of the n-dimensional
    cross-polytope onto a random hyperplane, for an odd dimension n.
    """
    try:
        n_str = input("Enter an odd integer dimension n (>=1): ")
        n = int(n_str)
        if n < 1 or n % 2 == 0:
            print("Error: The dimension n must be a positive odd integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # k = (n-1)/2
    k = (n - 1) // 2

    # Numerator of the formula: n*sqrt(n)
    numerator = n * math.sqrt(n)

    # Denominator of the formula: k! * pi^k
    k_factorial = math.factorial(k)
    pi_k = math.pow(math.pi, k)
    denominator = k_factorial * pi_k

    # Final result
    average_area = numerator / denominator

    # Output the components of the calculation as requested
    print(f"\nFor n = {n}:")
    print(f"The formula is n*sqrt(n) / (k! * pi^k), where k = (n-1)/2.")
    print(f"k = {k}")
    print(f"Numerator: n*sqrt(n) = {n}*sqrt({n}) = {numerator:.4f}")
    print(f"Denominator: k!*pi^k = {k}!*pi^{k} = {denominator:.4f}")
    
    print("\nResult:")
    print(f"The average projected area is: {average_area}")

if __name__ == '__main__':
    calculate_average_projection_area()