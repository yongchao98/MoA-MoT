import math

def calculate_average_projection_area():
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope,
    where n is an odd integer provided by the user.
    """
    try:
        # Prompt the user to enter the odd dimension n.
        n_str = input("Enter an odd dimension n (e.g., 3, 5, 7, ...): ")
        n = int(n_str)

        # Validate the input.
        if n < 3 or n % 2 == 0:
            print("Error: The dimension n must be an odd integer greater than or equal to 3.")
            return

        # For n = 2k + 1, we first find k.
        k = (n - 1) // 2

        # The derived formula for the average projected area is: sqrt(n) / (k! * (k-1)!)
        # We calculate each part of the formula.
        
        # Numerator
        numerator = math.sqrt(n)
        
        # Denominator
        k_factorial = math.factorial(k)
        k_minus_1_factorial = math.factorial(k - 1)
        denominator = k_factorial * k_minus_1_factorial

        # Final Result
        average_area = numerator / denominator

        # Print the results and the formula with the specific values for n.
        print(f"\nFor the {n}-dimensional cross-polytope:")
        print(f"The dimension n is {n}, which means k = (n-1)/2 = {k}.")
        print(f"The formula for the average projected area is: A = sqrt(n) / (k! * (k-1)!)")
        print(f"Substituting the values, we get: A = sqrt({n}) / (factorial({k}) * factorial({k-1}))")
        print(f"Which calculates to: A = {numerator:.5f} / ({k_factorial} * {k_minus_1_factorial})")
        print(f"The numerical value for the average projected area is: {average_area}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    calculate_average_projection_area()
