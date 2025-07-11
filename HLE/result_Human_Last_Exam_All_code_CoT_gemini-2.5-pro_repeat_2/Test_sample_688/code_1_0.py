import math

def calculate_cn(n):
    """
    Calculates the prefactor c_n for a given integer n >= 2.

    The formula for the prefactor is c_n = -(n-1) / n!
    
    Args:
        n (int): The number of nodes in the diagram.
    
    Returns:
        A tuple (numerator, denominator) representing the fraction.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer greater than or equal to 2.")
    
    numerator = -(n - 1)
    denominator = math.factorial(n)
    
    return numerator, denominator

def main():
    """
    Main function to calculate and display the prefactor c_n for several values of n.
    """
    print("The system-independent prefactor c_n is given by the formula: c_n = -(n-1) / n!")
    print("This code calculates c_n for a given n.")
    
    try:
        # We can calculate c_n for any integer n >= 2.
        # Let's demonstrate for n = 2, 3, 4, 5.
        for n in range(2, 6):
            num, den = calculate_cn(n)
            
            # As requested, outputting each number in the final equation.
            print(f"\nFor n = {n}:")
            print(f"c_{n} = -({n} - 1) / {n}!")
            
            # Find the simplified fraction
            common_divisor = math.gcd(abs(num), den)
            s_num = num // common_divisor
            s_den = den // common_divisor
            
            if s_den == 1:
                print(f"c_{n} = {s_num}")
            else:
                print(f"c_{n} = {s_num} / {s_den}")

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()