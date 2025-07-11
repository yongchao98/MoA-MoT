import math

def display_lower_bound():
    """
    This function prints the derived lower bound for the watermarking detection score E[S],
    and highlights the numbers used in the final formula as requested.
    """
    
    # Define the components of the formula as strings for display
    n = "n"
    alpha = "alpha"
    pi = "pi"
    
    # Define the integer numbers in the formula
    number_one = 1
    power_two = 2
    denominator_six = 6
    
    # Print the final expression for the lower bound
    print("A lower bound for the expected detection score E[S] is given by the expression:")
    print(f"E[S] > {n} * {alpha} + {n} * ({number_one} - {pi}^{power_two} / {denominator_six})")
    
    print("\nWhere:")
    print(f"- '{n}' is the number of tokens in the document.")
    print(f"- '{alpha}' is the average entropy of the model's output distribution per token.")
    print(f"- '{pi}' is the mathematical constant (approximately {math.pi}).")
    
    # As instructed, output each number present in the final equation
    print("\nThe specific numbers appearing in the final equation are:")
    print(f"1. The integer '{number_one}' in the term '({number_one} - ...)'.")
    print(f"2. The power '{power_two}' in the term '{pi}^{power_two}'.")
    print(f"3. The denominator '{denominator_six}' in the term '{pi}^{power_two} / {denominator_six}'.")

if __name__ == "__main__":
    display_lower_bound()