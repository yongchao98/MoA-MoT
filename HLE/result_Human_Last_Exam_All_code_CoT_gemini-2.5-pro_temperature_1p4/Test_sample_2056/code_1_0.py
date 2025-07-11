import numpy as np

def calculate_and_print_l_k(n: int):
    """
    Calculates the value of l_k(n) based on the derived analytical formula.

    Args:
        n (int): The dimension, must be >= 3.
    """
    if not isinstance(n, int) or n < 3:
        raise ValueError("n must be an integer greater than or equal to 3.")

    # Calculate the coefficients from the derived formula for the given n
    c1 = n - 1
    c2_num = 2 * n - 1
    c2_den = n
    c3 = 0.5
    c4 = n + 1
    
    # The constant k
    k = np.log(np.sqrt(2) + 1) # This is arccosh(sqrt(2))

    # Print the equation with numbers plugged in
    print(f"For n={n}, the exact expression for l_k(n) is:")
    print(f"l_k({n}) = ({c1})*ln(k) - k^2*({c2_num}/{c2_den}) + ({c3})*ln({c4})")
    print("\nwhere k = ln(sqrt(2) + 1).")
    
    # Calculate the final numerical value
    term1 = c1 * np.log(k)
    term2 = -(k**2) * (c2_num / c2_den)
    term3 = c3 * np.log(c4)
    result = term1 + term2 + term3
    
    # Print the numerical values
    print("\nThe numerical value is:")
    print(f"k = {k}")
    print(f"l_k({n}) = {result}")

if __name__ == '__main__':
    # Example usage for n=3
    calculate_and_print_l_k(n=3)
    print("\n" + "-"*30 + "\n")
    # Example usage for n=10
    calculate_and_print_l_k(n=10)
