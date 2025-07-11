def solve_bansenshukai_pattern():
    """
    This function decodes the Kunoichi no Jutsu pattern from the Bansenshukai.
    It interprets the pattern of black and white circles as a binary number,
    converts it to decimal, and prints the step-by-step calculation.
    """
    # The pattern is ⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤
    # Representing ⬤ as 1 and ○ as 0
    binary_string = "1001111011111"
    
    print(f"The pattern ⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤ can be represented in binary as: {binary_string}")
    
    decimal_value = 0
    equation_parts = []
    
    # Reverse the string to calculate powers of 2 from right to left (0 to n-1)
    reversed_binary = binary_string[::-1]
    
    for i in range(len(reversed_binary)):
        if reversed_binary[i] == '1':
            # Calculate the value for this position (1 * 2^i)
            part_value = 2**i
            decimal_value += part_value
            equation_parts.append(str(part_value))
            
    # Reverse the list of parts to match the original order of the binary string
    equation_parts.reverse()
    
    # Construct the final equation string
    equation_str = " + ".join(equation_parts)
    
    print("\nThe conversion to a decimal number is calculated as follows:")
    # Output each number in the final equation
    print(f"{equation_str} = {decimal_value}")

# Execute the function
solve_bansenshukai_pattern()