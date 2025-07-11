import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors and their exponents."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def multiplicative_order(a, n):
    """
    Calculates the multiplicative order of a modulo n.
    Assumes gcd(a, n) == 1.
    """
    if math.gcd(a, n) != 1:
        return -1 # Not defined
    
    result = 1
    value = a % n
    while value != 1:
        value = (value * a) % n
        result += 1
    return result

def is_np_possible_for_solvable_group(np, p):
    """
    Checks if np can be the number of Sylow p-subgroups in a solvable group
    based on Pazderski's theorem.
    """
    if np == 1:
        return True
    
    factors = get_prime_factorization(np)
    
    print(f"Checking y = {np}:")
    print(f"  Prime factorization of {np} is {factors}")
    
    for q, exponent in factors.items():
        if q == p:
            # The theorem applies to prime factors q != p.
            # n_p is not divisible by p.
            print(f"  Error: n_p ({np}) should not be divisible by p ({p}). This indicates an invalid n_p.")
            return False

        order_p_mod_q = multiplicative_order(p, q)
        
        print(f"  For prime factor q={q}:")
        print(f"    Exponent e = {exponent}")
        print(f"    Multiplicative order of p={p} modulo q={q} is {order_p_mod_q}")
        
        if exponent % order_p_mod_q != 0:
            print(f"    Condition failed: exponent {exponent} is not a multiple of order {order_p_mod_q}.")
            return False
        else:
            print(f"    Condition met: exponent {exponent} is a multiple of order {order_p_mod_q}.")
            
    return True

def find_minimum_y():
    """
    Finds the minimum y (number of Sylow 5-subgroups) that forces nonsolvability.
    """
    p = 5
    # y must be > 1 and y = 1 (mod 5)
    y = 6 
    while True:
        if not is_possible_np_for_solvable_group(y, p):
            print(f"\nFound the minimum value y = {y}.")
            print(f"A group with n_5 = {y} cannot be solvable, which guarantees it is nonsolvable.")
            print("This holds true regardless of the number of Sylow 3-subgroups.")
            return y
        y += 5

if __name__ == "__main__":
    min_y = find_minimum_y()
    print(f"\nFinal Answer: The minimum value of y is {min_y}.")
