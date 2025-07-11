import math

def calculate_cn(n):
    """
    Calculates the system-independent prefactor c_n for a given n.

    The formula for the prefactor is: c_n = -1 / (n * (n-2)!)

    This function prints the step-by-step calculation for clarity.
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer greater than or equal to 2. Got: {n}")
        return

    print(f"Calculating c_n for n = {n}:")
    
    # The final equation is c_n = -1 / (n * (n-2)!)
    # We will output each number/term involved in this equation.
    
    # Term 1: n
    print(f"  The term n is {n}")

    # Term 2: (n-2)!
    n_minus_2 = n - 2
    # The factorial() function requires a non-negative integer.
    if n_minus_2 < 0:
        print(f"  (n-2) is negative, factorial is not defined in the standard way.")
        return
    n_minus_2_factorial = math.factorial(n_minus_2)
    print(f"  The term (n-2)! is ({n}-2)! = {n_minus_2}! = {n_minus_2_factorial}")

    # Denominator: n * (n-2)!
    denominator = n * n_minus_2_factorial
    print(f"  The denominator n * (n-2)! is {n} * {n_minus_2_factorial} = {denominator}")

    # Final result: c_n
    cn = -1.0 / denominator
    cn_fraction = f"-1/{denominator}"
    print(f"  The prefactor c_n = {cn_fraction} = {cn}")
    print("-" * 25)

if __name__ == "__main__":
    # Calculate c_n for a few example values of n (from 2 to 5)
    for n_val in range(2, 6):
        calculate_cn(n_val)
