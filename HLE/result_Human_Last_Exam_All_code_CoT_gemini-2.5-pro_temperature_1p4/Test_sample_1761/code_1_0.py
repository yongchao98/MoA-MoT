import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1(2) on n-dimensional complex projective space.
    """
    try:
        n_str = input("Enter the dimension n of the complex projective space P^n: ")
        n = int(n_str)
        if n < 0:
            print("Error: The dimension n must be a non-negative integer.")
            return

        # The dimension is given by the formula C(n+1, 2) = n*(n+1)/2.
        
        # Calculate the components of the formula
        n_plus_1 = n + 1
        numerator = n * n_plus_1
        # Use integer division
        result = numerator // 2

        # Print the derivation with the specific numbers
        print(f"\nThe complex dimension is given by the formula h^0 = C(n+1, 2) = (n * (n + 1)) / 2")
        print(f"For n = {n}, the calculation is:")
        print(f"h^0 = ({n} * ({n} + 1)) / 2")
        print(f"h^0 = ({n} * {n_plus_1}) / 2")
        print(f"h^0 = {numerator} / 2")
        print(f"h^0 = {result}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer for n.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    calculate_dimension()