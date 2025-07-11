import math

def calculate_average_projection_area():
    """
    This script calculates the average area of a projection of the n-dimensional
    cross-polytope P onto a fixed hyperplane, after a uniformly random
    rotation of P. The dimension n is odd, n = 2k + 1.

    The average projection area is given by the formula:

      Area = sqrt(2k + 1) / k!

    This script will compute this value for a user-provided k.
    """
    try:
        k_str = input("Enter the non-negative integer k, where dimension n = 2k + 1: ")
        k = int(k_str)
        if k < 0:
            print("\nError: k must be a non-negative integer.")
            return

        n = 2 * k + 1

        numerator_val = math.sqrt(n)
        denominator_val = math.factorial(k)

        result = numerator_val / denominator_val

        print(f"\nFor the integer k = {k}, the odd dimension is n = 2*k + 1 = {n}.")
        print("The formula for the average projection area is sqrt(n) / k!.")
        print(f"We will compute sqrt({n}) / {k}!.")
        print("\n--- Calculation Steps ---")
        print(f"Numerator: sqrt({n}) = {numerator_val}")
        print(f"Denominator: {k}! = {denominator_val}")
        print("-----------------------")
        print(f"Final Result: Average Area = {numerator_val} / {denominator_val} = {result}")

    except ValueError:
        print("\nError: Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

calculate_average_projection_area()