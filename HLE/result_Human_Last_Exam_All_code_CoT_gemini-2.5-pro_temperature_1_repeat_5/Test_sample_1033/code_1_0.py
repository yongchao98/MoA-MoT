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

def find_next_prime_codes(start_code, count):
    """
    Finds the next 'count' codes in the sequence starting after 'start_code'.
    A code L1L2L3 is valid if the number v1*26^2 + v2*26 + v3 is prime.
    """
    v1 = ord(start_code[0]) - ord('A')
    v2 = ord(start_code[1]) - ord('A')
    v3 = ord(start_code[2]) - ord('A')

    found_codes = []
    
    print("Finding the next three codes in the sequence...\n")

    while len(found_codes) < count:
        # Increment the letter positions
        v3 += 1
        if v3 >= 26:
            v3 = 0
            v2 += 1
        if v2 >= 26:
            v2 = 0
            v1 += 1
        if v1 >= 26:
            # Reached the end of all possible codes
            break

        # Calculate the base-26 number
        num = v1 * 26**2 + v2 * 26 + v3

        if is_prime(num):
            code = chr(v1 + ord('A')) + chr(v2 + ord('A')) + chr(v3 + ord('A'))
            found_codes.append(code)
            print(f"Found code {len(found_codes)}: {code}")
            # Output the equation as requested
            print(f"Calculation: {v1} * 26^2 + {v2} * 26 + {v3} = {num}, which is a prime number.\n")

    return found_codes

if __name__ == "__main__":
    # The last code given in the sequence
    last_code = "NZX"
    
    # Find the next 3 codes
    next_three_codes = find_next_prime_codes(last_code, 3)
    
    # Format the final answer as requested
    final_answer = " ".join(next_three_codes)
    print(f"<<<{final_answer}>>>")