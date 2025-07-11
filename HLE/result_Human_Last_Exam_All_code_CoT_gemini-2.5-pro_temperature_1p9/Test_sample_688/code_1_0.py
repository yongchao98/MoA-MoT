import math

def calculate_cn_prefactor():
    """
    Asks the user for the number of particles 'n' and calculates the 
    prefactor c_n for the fully f-connected diagram in the virial expansion.
    
    The prefactor c_n is part of the n-th virial coefficient B_n, and its
    formula is derived as c_n = -(n-1) / n!
    """
    try:
        # Prompt user for input
        n_str = input("Enter the number of particles (n, an integer >= 2): ")
        n = int(n_str)

        # Validate the input
        if n < 2:
            print("Error: n must be an integer greater than or equal to 2.")
            return

        # Calculate the components of the formula
        numerator = -(n - 1)
        denominator = math.factorial(n)
        
        # Calculate the final result
        result = numerator / denominator

        # Print the derivation for the specific n
        print("\nThe prefactor c_n is given by the general formula: c_n = -(n-1) / n!")
        print(f"\nFor n = {n}, the calculation is:")
        
        # Output each number in the final equation
        print(f"c_{n} = -({n} - 1) / {n}!")
        print(f"c_{n} = {numerator} / {denominator}")
        print(f"c_{n} = {result}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_cn_prefactor()
