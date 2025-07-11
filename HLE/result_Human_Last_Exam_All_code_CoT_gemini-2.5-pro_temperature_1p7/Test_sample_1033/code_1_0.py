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
    """Gets the k-th prime number (1-based index)."""
    if k <= 0:
        return -1
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
    """
    # The next group starts with 'O', which is the 15th letter.
    n_char = 'O'
    n = ord(n_char) - ord('A') + 1
    
    # We need the next three terms, so k will be 1, 2, and 3.
    results = []
    print(f"Finding the next three terms for the group starting with '{n_char}' (n={n}).\n")
    
    for k in range(1, 4):
        p_k = get_kth_prime(k)
        
        # Calculate the value for the second letter
        v2_val = n + p_k
        v2_mod = (v2_val - 1) % 26 + 1
        char2 = chr(v2_mod - 1 + ord('A'))
        
        # Calculate the value for the third letter
        v3_val = n * p_k
        v3_mod = (v3_val - 1) % 26 + 1
        char3 = chr(v3_mod - 1 + ord('A'))
        
        result_str = f"{n_char}{char2}{char3}"
        results.append(result_str)
        
        print(f"Term k={k}:")
        print(f"  The {k}-st/nd/rd prime (p_{k}) is {p_k}.")
        print(f"  Second letter value: {n} + {p_k} = {v2_val}. Modulo 26, this is {v2_mod} ('{char2}').")
        print(f"  Third letter value: {n} * {p_k} = {v3_val}. Modulo 26, this is {v3_mod} ('{char3}').")
        print(f"  Result: {result_str}\n")
        
    print("The next three capital letters in the sequence are: " + " ".join(results))

solve()
