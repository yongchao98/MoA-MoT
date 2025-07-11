import math

def calculate_prefactor():
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    Ree-Hoover diagram contribution to the n-th virial coefficient.
    The user is prompted to enter the number of nodes, n.
    """
    try:
        n_str = input("Enter the number of nodes, n (must be an integer >= 2): ")
        n = int(n_str)

        if n < 2:
            print("\nError: The virial series starts at n=2. Please enter an integer greater than or equal to 2.")
            return

        print(f"\nDetermining the prefactor c_n for n = {n}.")

        # The derivation leads to the formula c_n = -(n-1)/n!
        print("The derived formula for the prefactor is: c_n = -(n - 1) / n!")

        num_val = n - 1
        den_val = math.factorial(n)

        # To simplify the fraction, we find the greatest common divisor.
        common_divisor = math.gcd(num_val, den_val)
        simplified_num = num_val // common_divisor
        simplified_den = den_val // common_divisor
        
        print(f"\nFor your value of n = {n}:")
        print(f"The numerator term (n - 1) is: {n} - 1 = {num_val}")
        print(f"The denominator term n! is: {den_val}")
        
        print(f"\nSubstituting these numbers into the formula:")
        # Show the equation with the computed values
        print(f"c_{n} = -({num_val}) / {den_val}")
        
        # Print simplified fraction
        if simplified_den == 1:
            print(f"The prefactor c_{n} is {-simplified_num}")
        else:
            # Show the simplified fraction
            print(f"As a simplified fraction, c_{n} = {-simplified_num} / {simplified_den}")
            
        # Print floating point result
        result = -num_val / den_val
        print(f"As a decimal, c_{n} = {result}")

    except ValueError:
        print("\nError: Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"\nAn error occurred: {e}")

calculate_prefactor()