import sys

def get_shapley_formula(n):
    """
    Calculates and prints the formula for the Shapley value c_k for player p_k.

    Args:
        n (int): The total number of people in the band.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.")
        return

    # Calculate the sum of the first n integers (S1)
    s1 = n * (n + 1) // 2

    # Calculate the sum of the first n squares (S2)
    s2 = n * (n + 1) * (2 * n + 1) // 6

    # The formula for c_k is a quadratic in k: A*k^2 + B*k
    # c_k = -s1^2 * k^2 + (s1^3 + s1*s2) * k
    
    # Coefficient for k^2
    coeff_k2 = -s1**2
    
    # Coefficient for k
    coeff_k1 = s1**3 + s1 * s2

    print(f"For n = {n}, the fair amount of money c_k for player p_k is given by the formula:")
    # The final formula is c_k = (-S1^2)k^2 + (S1^3 + S1*S2)k
    # We output each number in the equation.
    print(f"c_k = {coeff_k2}*k^2 + {coeff_k1}*k")

if __name__ == '__main__':
    # You can change the value of n here to see the formula for different group sizes.
    # Example usage:
    try:
        if len(sys.argv) > 1:
            n_value = int(sys.argv[1])
        else:
            # Default value if no command-line argument is given.
            n_value = 4
            print(f"No value for n provided. Using default n = {n_value}.")
            
        get_shapley_formula(n_value)

    except (ValueError, IndexError):
        print("Invalid input. Please provide an integer value for n.")
        print("Usage: python your_script_name.py <n>")