import collections

def multiply(g1, g2):
    """
    Multiplies two elements in the generalized quaternion group Q128.
    Elements are represented as tuples (type, exponent), where:
    - type 0: element is of the form x^exponent
    - type 1: element is of the form y*x^exponent
    The group relations are x^64 = 1, y^2 = x^32, y*x*y^-1 = x^-1.
    """
    t1, e1 = g1
    t2, e2 = g2
    N = 64  # Order of the cyclic subgroup <x>
    Z = 32  # Power of x for y^2

    if t1 == 0 and t2 == 0:
        # x^e1 * x^e2 = x^(e1+e2)
        return (0, (e1 + e2) % N)
    elif t1 == 0 and t2 == 1:
        # x^e1 * (y*x^e2) = y*x^(-e1)*x^e2 = y*x^(e2-e1)
        return (1, (e2 - e1) % N)
    elif t1 == 1 and t2 == 0:
        # (y*x^e1) * x^e2 = y*x^(e1+e2)
        return (1, (e1 + e2) % N)
    elif t1 == 1 and t2 == 1:
        # (y*x^e1) * (y*x^e2) = y*x^e1*y*x^e2 = y*y*x^(-e1)*x^e2 = y^2*x^(e2-e1) = x^32*x^(e2-e1)
        return (0, (Z + e2 - e1) % N)

def power(g, k):
    """
    Computes g^k for an element g in Q128 using exponentiation by squaring.
    """
    if k < 0:
        raise ValueError("Power must be a non-negative integer.")
    
    # Identity element is x^0
    identity = (0, 0)
    
    if k == 0:
        return identity
    
    res = identity
    base = g
    while k > 0:
        if k % 2 == 1:
            res = multiply(res, base)
        base = multiply(base, base)
        k //= 2
    return res

def solve():
    """
    Calculates the number of power subgroups in Q128.
    """
    N = 64
    Q128 = []
    for i in range(N):
        Q128.append((0, i))
    for i in range(N):
        Q128.append((1, i))

    odd_k_subgroups = set()
    even_k_subgroups = set()
    
    # The exponent of Q128 is 64. We check k from 1 to 64.
    # G^k = G^j if gcd(k, exp(G)) = gcd(j, exp(G)).
    for k in range(1, N + 1):
        # A frozenset is an immutable, hashable set.
        H_k = frozenset(power(g, k) for g in Q128)
        if k % 2 != 0:
            odd_k_subgroups.add(H_k)
        else:
            even_k_subgroups.add(H_k)
            
    num_odd = len(odd_k_subgroups)
    num_even = len(even_k_subgroups)
    
    # The set of subgroups from odd k is disjoint from the set from even k,
    # as G is non-cyclic, while subgroups from even k are all cyclic.
    total = num_odd + num_even

    print("The generalized quaternion group of size 128 is Q128.")
    print("A power subgroup G^k is the set {g^k | g in G}.")
    print("We compute these subgroups for exponents k from 1 up to 64 (the exponent of Q128).")
    print(f"For any odd integer k, G^k is the group G itself. This gives {num_odd} unique subgroup.")
    print(f"For even integers k (from 2 to 64), we find {num_even} distinct subgroups.")
    print("The group G itself is not found among the subgroups for even k.")
    print(f"Total number of power subgroups = {num_odd} + {num_even} = {total}")

solve()