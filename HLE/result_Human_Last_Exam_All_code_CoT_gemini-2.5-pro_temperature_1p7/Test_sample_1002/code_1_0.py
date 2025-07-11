import sys

def calculate_limit_exponent(k):
    """
    Calculates the limit lim_{m->inf} ln(f(m))/ln(m) for a given integer k >= 2.
    
    The limit is given by the formula (2k - 2) / (2k - 1).
    This function demonstrates the calculation and prints the components of the final equation.
    
    Args:
        k (int): An integer greater than or equal to 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: Input 'k' must be an integer greater than or equal to 2.")
        return

    # Define the numbers in the numerator's equation
    num_coeff = 2
    num_sub = 2

    # Define the numbers in the denominator's equation
    den_coeff = 2
    den_sub = 1

    # Calculate numerator and denominator
    numerator = num_coeff * k - num_sub
    denominator = den_coeff * k - den_sub

    # Calculate the final result
    result = numerator / denominator

    # Output the explanation and each number in the final equation
    print(f"For k = {k}, we are calculating the value of the limit, which is given by the formula (2*k - 2)/(2*k - 1).")
    print("\n--- Calculation Steps ---")
    
    # Printing the numerator calculation showing all numbers
    print("Numerator calculation:")
    print(f"  ({num_coeff} * {k}) - {num_sub} = {numerator}")

    # Printing the denominator calculation showing all numbers
    print("Denominator calculation:")
    print(f"  ({den_coeff} * {k}) - {den_sub} = {denominator}")
    
    print("\n--- Final Equation ---")
    # Printing the final equation with the calculated numbers
    print(f"  {numerator} / {denominator} = {result}")

if __name__ == '__main__':
    # You can change the value of k here.
    # It must be an integer >= 2.
    # Let's use k=4 as an example.
    k_value = 4
    try:
        # Check if a command-line argument is provided
        if len(sys.argv) > 1:
            k_value = int(sys.argv[1])
    except (ValueError, IndexError):
        print(f"Invalid input. Using default k = {k_value}.")
        
    calculate_limit_exponent(k_value)
