import roman

def get_roman_ascii_sum(n):
    """
    Calculates the sum of ASCII values for the characters
    in the Roman numeral representation of a number.
    """
    roman_numeral = roman.toRoman(n)
    
    ascii_sum = 0
    equation_parts = []
    
    for char in roman_numeral:
        ascii_val = ord(char)
        ascii_sum += ascii_val
        equation_parts.append(f"ord('{char}')")
        
    return ascii_sum, roman_numeral, equation_parts

def solve_sequence():
    """
    Finds the next number in the sequence.
    """
    # The given sequence has 33 numbers, so we need to find the 34th.
    n = 34
    
    # Calculate the sum and get the parts for the equation
    next_number, roman_n, equation = get_roman_ascii_sum(n)
    
    # As requested, print each number in the final equation
    equation_string = " + ".join(equation)
    value_string = " + ".join([str(ord(c)) for c in roman_n])
    
    print(f"To find the next number in the sequence, we calculate it for n = {n}.")
    print(f"The Roman numeral for {n} is '{roman_n}'.")
    print("The value is the sum of the ASCII values of its characters.")
    print(f"\nFinal Equation:")
    print(f"{equation_string} = {value_string} = {next_number}")

solve_sequence()
<<<423>>>