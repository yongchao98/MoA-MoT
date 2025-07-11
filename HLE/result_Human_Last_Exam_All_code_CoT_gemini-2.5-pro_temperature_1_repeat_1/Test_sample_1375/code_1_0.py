import math

def calculate_average_projection_area(n):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope.
    The dimension n must be an odd integer >= 3.
    """
    if not isinstance(n, int) or n < 3 or n % 2 == 0:
        print("Error: Input must be an odd integer greater than or equal to 3.")
        return

    # For n = 2k + 1, k = (n-1)/2
    k = (n - 1) // 2

    # The formula for the average area is sqrt(n) / (k! * (k-1)!)
    try:
        k_factorial = math.factorial(k)
        k_minus_1_factorial = math.factorial(k - 1)
    except ValueError:
        print("Error: Factorial calculation failed. k must be non-negative.")
        return

    numerator = math.sqrt(n)
    denominator = k_factorial * k_minus_1_factorial
    average_area = numerator / denominator

    # Output the details of the calculation
    print(f"\nFor the {n}-dimensional cross-polytope:")
    print(f"The dimension is n = {n}, which is odd.")
    print(f"We use the formula: A_avg = sqrt(n) / (k! * (k-1)!), where k = (n-1)/2.")
    print(f"For n = {n}, k = ( {n} - 1 ) / 2 = {k}.")
    print("\nCalculation steps:")
    print(f"Numerator = sqrt({n}) = {numerator}")
    print(f"Denominator = {k}! * ({k}-1)! = {k_factorial} * {k_minus_1_factorial} = {denominator}")
    print(f"Average Area = {numerator} / {denominator}")
    print(f"Result: {average_area}")

def main():
    """
    Main function to get user input and call the calculation function.
    """
    try:
        n_input = input("Enter the odd dimension n (must be >= 3): ")
        n = int(n_input)
        calculate_average_projection_area(n)
    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
