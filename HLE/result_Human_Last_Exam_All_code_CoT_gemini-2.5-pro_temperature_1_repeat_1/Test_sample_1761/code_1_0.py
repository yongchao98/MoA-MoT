import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1(2) on complex projective space P^n.
    """
    try:
        n_str = input("Enter the dimension 'n' of the complex projective space P^n (a positive integer): ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension 'n' must be a positive integer.")
            return

        # The dimension is given by the formula n * (n + 1) / 2.
        # This is the dimension of the space of (n+1)x(n+1) skew-symmetric matrices.
        numerator = n * (n + 1)
        dimension = numerator // 2
        
        # Print the explanation and the final calculation step-by-step.
        print(f"\nFor P^{n} with n = {n}:")
        print("The complex dimension of H^0(P^n, Omega^1(2)) is given by the formula: (n * (n + 1)) / 2")
        print("This corresponds to the number of ways to choose 2 items from n+1, C(n+1, 2).")
        
        print("\nCalculation:")
        print(f"({n} * ({n} + 1)) / 2 = ({n} * {n+1}) / 2 = {numerator} / 2 = {dimension}")
        
        print(f"\nThe complex dimension is {dimension}.")

    except ValueError:
        print("Error: Invalid input. Please enter a valid integer for 'n'.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_dimension()