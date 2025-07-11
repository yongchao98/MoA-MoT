def solve_bansenshukai_riddle():
    """
    Analyzes the Kunoichi no Jutsu circle pattern from the Bansenshukai
    by interpreting it as a binary number and printing the conversion equation.
    
    The pattern is: ⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤
    Let ⬤ = 1 and ○ = 0.
    """
    
    # The pattern from the scroll
    circle_pattern = "⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤"
    
    # Convert the pattern to a binary string
    binary_string = circle_pattern.replace('⬤', '1').replace('○', '0')
    
    # Initialize variables for the equation
    decimal_value = 0
    equation_parts = []
    num_bits = len(binary_string)
    
    # Calculate the decimal value and build the equation string
    # Iterate in reverse to match powers of 2 (2^0, 2^1, etc.)
    for i, bit in enumerate(reversed(binary_string)):
        power = i
        if bit == '1':
            term_value = 2**power
            decimal_value += term_value
            equation_parts.append(str(term_value))
        else:
            # Add 0 for the empty circles to show all parts of the equation
            equation_parts.append("0")
            
    # The equation should read from highest power to lowest
    equation_parts.reverse()
    
    equation_string = " + ".join(equation_parts)
    
    print(f"Interpreting the pattern as a binary number (⬤=1, ○=0):")
    print(f"Binary: {binary_string}")
    print("\nThis corresponds to the following equation:")
    # We print each number in the final equation as requested
    print(f"{equation_string} = {decimal_value}")

solve_bansenshukai_riddle()