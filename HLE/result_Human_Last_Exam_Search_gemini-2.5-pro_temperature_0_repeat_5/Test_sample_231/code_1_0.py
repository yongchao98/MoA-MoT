import sympy

def find_prime_relation(limit_q):
    for q in sympy.primerange(3, limit_q):
        p_q = sympy.prime(q)
        if (p_q + 1) % q == 0:
            k = (p_q + 1) / q
            if k > 2:
                n = p_q
                print(f"Found a solution: q = {q}, p_q = {n}, k = {k}")
                return k, n
    return None, None

# Let's run this with a reasonable limit.
# find_prime_relation(1000)