import math

def is_prime(n):
    """
    Checks if a number is prime using an efficient algorithm.
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

def num_to_code(n):
    """
    Converts a decimal number to its 3-letter base-26 representation.
    """
    if not (0 <= n < 26**3):
        return None  # Out of range for a 3-letter code
    l1 = n // (26**2)
    rem = n % (26**2)
    l2 = rem // 26
    l3 = rem % 26
    return chr(ord('A') + l1) + chr(ord('A') + l2) + chr(ord('A') + l3)

def code_to_num(code):
    """
    Converts a 3-letter base-26 code to its decimal number.
    """
    l1 = ord(code[0]) - ord('A')
    l2 = ord(code[1]) - ord('A')
    l3 = ord(code[2]) - ord('A')
    return l1 * (26**2) + l2 * 26 + l3

def find_next_codes(start_code, count=3):
    """
    Finds the next 'count' prime-based codes after the start_code.
    """
    start_num = code_to_num(start_code)
    found_codes = []
    current_num = start_num + 1

    while len(found_codes) < count:
        if is_prime(current_num):
            code = num_to_code(current_num)
            if code:
                found_codes.append(code)
        current_num += 1
    
    return found_codes

if __name__ == "__main__":
    # The last code given in the sequence
    last_code = 'NZX'
    
    # Find the next three codes
    next_three_codes = find_next_codes(last_code, 3)
    
    print(f"The sequence is based on prime numbers represented in base-26 (A=0 to Z=25).")
    print(f"The last code provided is {last_code}, which corresponds to the prime number {code_to_num(last_code)}.")
    print("The next three codes in the sequence are:")
    print("-" * 30)

    # Print the results with the requested equation format
    for code in next_three_codes:
        num = code_to_num(code)
        l1 = ord(code[0]) - ord('A')
        l2 = ord(code[1]) - ord('A')
        l3 = ord(code[2]) - ord('A')
        print(f"Code: {code}")
        # Printing the equation with each number as requested
        print(f"Equation: {l1} * 26^2 + {l2} * 26 + {l3} = {num}")
        print("-" * 30)