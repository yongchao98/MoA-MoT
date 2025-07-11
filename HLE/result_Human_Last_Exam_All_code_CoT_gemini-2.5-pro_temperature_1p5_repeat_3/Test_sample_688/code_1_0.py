import math

def calculate_prefactor_cn(n):
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    contribution to the n-th virial coefficient.

    The formula is c_n = -(n-1) / n!

    Args:
        n (int): The number of particles (must be >= 2).

    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer and greater than or equal to 2.")
    
    numerator = n - 1
    denominator = math.factorial(n)
    
    return -numerator / denominator

if __name__ == "__main__":
    print("The prefactor c_n is given by the formula: c_n = -(n - 1) / n!")
    print("Below are the values for n from 2 to 10.")
    print("-" * 40)

    for n_val in range(2, 11):
        # Numbers in the equation
        n_numerator = n_val - 1
        n_denominator = math.factorial(n_val)
        result = -n_numerator / n_denominator
        
        # Outputting each number in the final equation as requested
        print(f"For n = {n_val}:")
        print(f"c_{n_val} = -({n_val} - 1) / {n_val}!")
        print(f"   = -{n_numerator} / {n_denominator}")
        print(f"   = {result}")
        print("-" * 40)