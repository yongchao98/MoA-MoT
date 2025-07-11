import math

def is_prime(n):
    """
    Checks if a number is prime using an optimized trial division method.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def code_to_num(code):
    """
    Converts a 3-letter code (as base 26) to a number (base 10).
    A=0, B=1, ..., Z=25.
    """
    if len(code) != 3:
        raise ValueError("Code must be three letters long.")
    n1 = ord(code[0]) - ord('A')
    n2 = ord(code[1]) - ord('A')
    n3 = ord(code[2]) - ord('A')
    return n1 * 26**2 + n2 * 26 + n3

def num_to_code(n):
    """
    Converts a number (base 10) to its 3-letter code representation (base 26).
    Also returns the component values for the equation.
    """
    val1 = n // (26**2)
    remainder = n % (26**2)
    val2 = remainder // 26
    val3 = remainder % 26
    
    char1 = chr(val1 + ord('A'))
    char2 = chr(val2 + ord('A'))
    char3 = chr(val3 + ord('A'))
    
    code = f"{char1}{char2}{char3}"
    return code, val1, val2, val3

def find_next_codes_in_sequence():
    """
    Finds and prints the next three codes in the prime number sequence.
    """
    # The last code provided in the sequence
    last_code = "NZX"
    
    # Convert the last code to its numeric value to know where to start searching
    start_num = code_to_num(last_code)
    
    print(f"The sequence consists of 3-letter representations of consecutive prime numbers in base 26.")
    print(f"The last term is {last_code}, which corresponds to the number {start_num}.")
    print("Finding the next three prime numbers...\n")
    
    found_codes = []
    num_to_check = start_num + 1
    
    while len(found_codes) < 3:
        if is_prime(num_to_check):
            code, v1, v2, v3 = num_to_code(num_to_check)
            # Format the output to show the code and the equation that generates it
            equation = f"{code}: {v1} * 26^2 + {v2} * 26 + {v3} = {num_to_check}"
            found_codes.append(equation)
        num_to_check += 1
        
    # Print the final results
    for result in found_codes:
        print(result)

# Execute the main function
find_next_codes_in_sequence()