import math

def is_prime(num):
    """Checks if a number is prime."""
    if num < 2:
        return False
    for i in range(2, int(math.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True

def get_kth_prime(k):
    """Finds the k-th prime number."""
    count = 0
    num = 1
    while count < k:
        num += 1
        if is_prime(num):
            count += 1
    return num

def solve():
    """
    Solves the sequence puzzle to find the next three terms.
    The rule is determined to be based on the group index (n),
    term index within the group (k), and the k-th prime number (p_k).
    """
    # The next group is 'O', which is the 15th letter of the alphabet.
    n = 15
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    next_terms = []

    print("The next three capital letters in the sequence are:")
    # We need to find the next 3 terms, so for k=1, k=2, and k=3.
    for k in range(1, 4):
        p_k = get_kth_prime(k)

        # The rule for the sequence:
        # L1 is the n-th letter.
        # L2_val = (p_k + n - 1) mod 26
        # L3_val = (p_k * n - 1) mod 26
        # Using A=0, B=1, ... Z=25 for values.
        # n is 1-based (A=1, ... O=15)
        
        l1_char = alphabet[n - 1]
        
        l2_val = (p_k + n - 1) % 26
        l2_char = alphabet[l2_val]
        
        l3_val = (p_k * n - 1) % 26
        l3_char = alphabet[l3_val]
        
        term = f"{l1_char}{l2_char}{l3_char}"
        next_terms.append(term)
    
    print(", ".join(next_terms))

solve()