import math

def is_prime(n):
    """Checks if a number is prime."""
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

def find_next_codes():
    """Finds the next three prime-representing codes after NZX."""
    # Start numerically from the code 'NZX' where N=13, Z=25, X=23.
    n1, n2, n3 = 13, 25, 23
    
    found_codes = []
    print("Finding the next three codes in the sequence after NZX...")
    print("-" * 60)
    print("The rule: A code L1L2L3 is in the sequence if its base-26 value is a prime number.")
    print("Value = (L1 position) * 26^2 + (L2 position) * 26 + (L3 position)")
    print("-" * 60)

    while len(found_codes) < 3:
        # Increment to the next code alphabetically
        n3 += 1
        if n3 > 25:
            n3 = 0
            n2 += 1
            if n2 > 25:
                n2 = 0
                n1 += 1
                if n1 > 25:
                    break # Should not happen

        # Calculate the base-26 value
        value = n1 * 26**2 + n2 * 26 + n3
        
        # Check if the value is prime
        if is_prime(value):
            code = chr(n1 + ord('A')) + chr(n2 + ord('A')) + chr(n3 + ord('A'))
            found_codes.append(code)
            # Output the equation as requested
            print(f"Found code: {code}")
            print(f"Equation: {n1} * 26^2 + {n2} * 26 + {n3} = {value} (which is a prime number)\n")

find_next_codes()
<<<NZZ OAD OAF>>>