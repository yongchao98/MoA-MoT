import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
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

def find_next_codes(last_code, count):
    """
    Finds the next 'count' codes in the sequence starting after 'last_code'.
    """
    # Convert the last code's letters to 0-indexed values
    v1 = ord(last_code[0]) - ord('A')
    v2 = ord(last_code[1]) - ord('A')
    v3 = ord(last_code[2]) - ord('A')

    found_codes = []
    
    # Iterate through all possible codes in lexicographical order
    while len(found_codes) < count:
        # Increment to the next code
        v3 += 1
        if v3 > 25:
            v3 = 0
            v2 += 1
            if v2 > 25:
                v2 = 0
                v1 += 1
        
        # This shouldn't happen for this problem, but as a safeguard
        if v1 > 25:
            break

        # Calculate the base-26 integer value
        num = v1 * 26**2 + v2 * 26 + v3
        
        # Check if the number is prime
        if is_prime(num):
            c1 = chr(ord('A') + v1)
            c2 = chr(ord('A') + v2)
            c3 = chr(ord('A') + v3)
            code = f"{c1}{c2}{c3}"
            
            equation = f"{c1}{c2}{c3} = {v1} * 26^2 + {v2} * 26 + {v3} = {v1 * 26**2} + {v2 * 26} + {v3} = {num}"
            found_codes.append((code, equation))

    return found_codes

def main():
    """
    Main function to solve the sequence puzzle.
    """
    # The last code from the provided sequence
    last_code_in_sequence = "NZX"
    
    # We need to find the next 3 codes
    num_to_find = 3
    
    print(f"The sequence is formed by 3-letter codes that represent prime numbers in base 26.")
    print(f"Searching for the next {num_to_find} codes after {last_code_in_sequence}...")
    print("-" * 20)
    
    next_codes = find_next_codes(last_code_in_sequence, num_to_find)
    
    final_answer = []
    for code, equation in next_codes:
        print(f"Found code: {equation}, which is prime.")
        final_answer.append(code)
    
    print("-" * 20)
    print(f"The next three capital letters in the sequence are: {', '.join(final_answer)}")


if __name__ == "__main__":
    main()
