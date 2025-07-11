def is_prime(n):
    """
    Checks if a number is prime using an efficient method.
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

def to_base26_triplet(p):
    """
    Converts a number into a base-26 triplet (n1, n2, n3).
    """
    n3 = p % 26
    p //= 26
    n2 = p % 26
    p //= 26
    n1 = p
    return n1, n2, n3

def to_letters(n1, n2, n3):
    """
    Converts a base-26 triplet to a three-letter string.
    """
    c1 = chr(ord('A') + n1)
    c2 = chr(ord('A') + n2)
    c3 = chr(ord('A') + n3)
    return f"{c1}{c2}{c3}"

def solve():
    """
    Finds the next three terms in the sequence.
    """
    # The last term given is NZX.
    # N=13, Z=25, X=23
    # In base 10, this is 13*26^2 + 25*26 + 23 = 9461.
    last_val = 9461
    
    # The next terms must start with 'O', which is 14 in 0-indexed alphabet.
    target_n1 = 14
    
    found_terms = []
    
    # We start searching for primes from the number after the last term's value.
    current_num = last_val + 1
    
    print("The next three terms are found by finding the next prime numbers in the sequence.")
    print("The sequence encodes prime numbers in base 26 (A=0, Z=25).")
    print("Last term NZX corresponds to 13*26^2 + 25*26 + 23 = 9461.")
    print("Next terms must start with 'O' (14), so we search for primes p > 9461 where p / 26^2 = 14.\n")

    while len(found_terms) < 3:
        if is_prime(current_num):
            n1, n2, n3 = to_base26_triplet(current_num)
            
            if n1 == target_n1:
                term_str = to_letters(n1, n2, n3)
                found_terms.append(term_str)
                # Output the equation for each found term
                print(f"Found term: {term_str}")
                print(f"Equation: {n1} * 26^2 + {n2} * 26 + {n3} = {current_num} (which is prime)")
                
        current_num += 1

    print("\nThe three capital letters that come next are:", ", ".join(found_terms))

solve()
<<<OAD, OAJ, OAP>>>