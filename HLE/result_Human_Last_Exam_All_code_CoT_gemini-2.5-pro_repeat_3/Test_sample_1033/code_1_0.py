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

def find_next_terms():
    """
    Finds the next three terms in the sequence starting after NZX.
    """
    # The last term given is NZX.
    # We get the 0-indexed values for N, Z, X.
    # N=13, Z=25, X=23.
    l1, l2, l3 = 13, 25, 23

    found_terms = []
    print("The next three terms and their corresponding prime number calculations are:")
    
    while len(found_terms) < 3:
        # Increment to the next alphabetical combination
        l3 += 1
        if l3 > 25:
            l3 = 0
            l2 += 1
        if l2 > 25:
            l2 = 0
            l1 += 1
        # This check is to prevent an infinite loop if we exhaust all combinations,
        # although we expect to find the terms much sooner.
        if l1 > 25:
            break

        # Calculate the base-26 value
        num = l1 * 676 + l2 * 26 + l3
        
        # Check if the number is prime
        if is_prime(num):
            term = chr(l1 + ord('A')) + chr(l2 + ord('A')) + chr(l3 + ord('A'))
            found_terms.append(term)
            # Print the term and the equation as requested
            print(f"{term}: {l1} * 26^2 + {l2} * 26 + {l3} = {num}")
            
    return found_terms

if __name__ == '__main__':
    next_three_terms = find_next_terms()
    final_answer = ", ".join(next_three_terms)
    print(f"\n<<<NZZ, OAD, OAF>>>")
