import math

def calculate_cn(n):
    """
    Calculates the system-independent prefactor c_n for the virial expansion.
    
    The formula for the prefactor is c_n = -(n-1) / n!
    
    Args:
        n (int): The number of particles in the cluster (must be >= 2).
        
    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer greater than or equal to 2.")
    
    numerator = -(n - 1)
    denominator = math.factorial(n)
    
    return float(numerator) / denominator

def main():
    """
    Main function to demonstrate the calculation of c_n for a range of n.
    """
    print("Calculating the prefactor c_n for the fully f-connected Ree-Hoover diagram.")
    print("The formula is: c_n = -(n-1) / n!")
    print("-" * 50)
    
    for n_val in range(2, 9):
        num = -(n_val - 1)
        den = math.factorial(n_val)
        result = calculate_cn(n_val)
        
        # Output each number in the final equation as requested.
        print(f"For n = {n_val}:")
        print(f"  c_{n_val} = -({n_val}-1) / {n_val}! = {num} / {den} = {result}")

if __name__ == "__main__":
    main()
