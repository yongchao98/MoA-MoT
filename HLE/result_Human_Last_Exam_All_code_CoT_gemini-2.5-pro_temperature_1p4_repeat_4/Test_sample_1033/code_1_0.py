import math

def is_prime(n):
    """
    Checks if a number is prime using trial division.
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

def code_to_number(code):
    """
    Converts a 3-letter code (e.g., 'ACB') to its base-10 value.
    A=0, B=1, ... Z=25.
    """
    c1 = ord(code[0]) - ord('A')
    c2 = ord(code[1]) - ord('A')
    c3 = ord(code[2]) - ord('A')
    return c1 * 26**2 + c2 * 26 + c3

def number_to_code(n):
    """
    Converts a base-10 number to its 3-letter base-26 code.
    """
    c1_val = n // (26**2)
    remainder = n % (26**2)
    c2_val = remainder // 26
    c3_val = remainder % 26
    
    c1 = chr(c1_val + ord('A'))
    c2 = chr(c2_val + ord('A'))
    c3 = chr(c3_val + ord('A'))
    
    return f"{c1}{c2}{c3}", c1_val, c2_val, c3_val

def find_next_terms():
    """
    Finds the next three terms in the sequence.
    """
    last_term = "NZX"
    current_number = code_to_number(last_term)
    
    found_terms = []
    
    num_to_find = 3
    candidate = current_number + 1
    
    while len(found_terms) < num_to_find:
        if is_prime(candidate):
            code, c1, c2, c3 = number_to_code(candidate)
            # Store the term and its components for final printout
            found_terms.append({
                "code": code,
                "prime": candidate,
                "c1": c1,
                "c2": c2,
                "c3": c3
            })
        candidate += 1
        
    # Print the final results
    print("The next three capital letters in the sequence are:")
    final_answer = []
    for term in found_terms:
        print(f"Term: {term['code']}")
        print(f"This corresponds to the prime number {term['prime']}.")
        print(f"Equation: {term['prime']} = {term['c1']} * 26^2 + {term['c2']} * 26 + {term['c3']}")
        print("-" * 20)
        final_answer.append(term['code'])
    
    return final_answer

if __name__ == '__main__':
    result = find_next_terms()
    # Required final answer format
    print(f"<<<{', '.join(result)}>>>")
