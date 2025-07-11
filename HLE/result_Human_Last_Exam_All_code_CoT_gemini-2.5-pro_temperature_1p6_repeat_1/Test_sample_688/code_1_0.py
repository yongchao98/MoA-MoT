import math

def calculate_and_print_cn(n):
    """
    Calculates and prints the prefactor c_n for a given n.

    The formula for c_n is derived as -(n-1) / n!
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer greater than or equal to 2, but got {n}")
        return

    # Calculate the components of the formula
    numerator = n - 1
    denominator = math.factorial(n)
    
    # Calculate c_n
    cn_value = -numerator / denominator

    # Print the results, showing each number in the equation
    print(f"For n = {n}:")
    print(f"  The formula for the prefactor is c_n = -(n - 1) / n!")
    print(f"  Plugging in n = {n}:")
    print(f"  c_{n} = -({numerator}) / ({denominator})")
    print(f"  c_{n} = {cn_value}")
    print("-" * 25)

# Demonstrate the calculation for n = 2, 3, 4, 5
if __name__ == "__main__":
    calculate_and_print_cn(2)
    calculate_and_print_cn(3)
    calculate_and_print_cn(4)
    calculate_and_print_cn(5)