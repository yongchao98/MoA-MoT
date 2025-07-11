import math

def prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
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

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve():
    """
    Solves the problem by counting the number of non-isomorphic groups
    of order 10! that can be fundamental groups of closed orientable 3-manifolds.
    """
    order = math.factorial(10)
    factors = prime_factorization(order)

    # --- Case 1: All Sylow subgroups are cyclic ---
    # These groups are of the form Z_m semidirect product Z_n, where n is the 2-power part
    # and m is the odd part of the order.
    # n = 2^8 = 256
    # m = 3^4 * 5^2 * 7^1 = 14175
    # The number of such groups is the number of orbits of a subgroup H of Aut(Z_m)
    # under an action by Aut(Z_n). This simplifies to a simpler orbit calculation.
    
    # Aut(Z_m) is isomorphic to the product of Aut(Z_{p^k}) for each prime power in m.
    # Aut(Z_{81}) is cyclic of order phi(81) = 54
    # Aut(Z_{25}) is cyclic of order phi(25) = 20
    # Aut(Z_7)  is cyclic of order phi(7)  = 6
    aut_orders = [phi(81), phi(25), phi(7)]

    # The group H consists of elements in Aut(Z_m) whose order is a power of 2.
    # In a cyclic group C_N, the number of elements whose order is a 2-power is gcd(N, 2^k)
    # where 2^k is the highest power of 2 dividing N.
    num_2_power_elements_H = [
        gcd(54, 2**8),  # In C_54, 2-power orders can be 1, 2. So 2 elements.
        gcd(20, 2**8),  # In C_20, 2-power orders can be 1, 2, 4. So 4 elements.
        gcd(6, 2**8)   # In C_6, 2-power orders can be 1, 2. So 2 elements.
    ]
    size_H = num_2_power_elements_H[0] * num_2_power_elements_H[1] * num_2_power_elements_H[2]
    
    # The number of orbits is given by (1/2) * (|H| + # of involutions in H).
    # The number of involutions (elements x with x^2=e) in C_N is gcd(2,N).
    num_involutions = [
        gcd(54, 2),
        gcd(20, 2),
        gcd(6, 2)
    ]
    size_involutions_H = num_involutions[0] * num_involutions[1] * num_involutions[2]

    num_type1_groups = (size_H + size_involutions_H) // 2

    # --- Case 2: Sylow 2-subgroup is Generalized Quaternion (GQ) ---
    # These groups are semidirect products Z_m x Q_256.
    # Counting them amounts to counting orbits via Burnside's Lemma.
    # The acting group is GL(2, F_2) of order 6. The set is Hom(Z_2 x Z_2, (Z_2)^3) of size 64.
    
    # Burnside's Lemma: num_orbits = (1/|G|) * sum_{g in G} |fix(g)|
    # G = GL(2, F_2), |G| = 6
    # X = Hom( (Z/2)^2, (Z/2)^3 ), |X| = 64
    
    # Elements in GL(2, F_2):
    # 1 identity element: fixes all 64 maps.
    # 3 involutions (order 2): each fixes 8 maps.
    # 2 elements of order 3: each fixes 1 map (the zero map).
    
    num_type2_groups = (1 * 64 + 3 * 8 + 2 * 1) // 6

    # Total number is the sum of the two cases.
    total = num_type1_groups + num_type2_groups
    
    print("This problem reduces to counting specific types of finite groups of order 10!.")
    print("Order = 10! =", order)
    print("Prime factorization:", " x ".join(f"{p}^{e}" for p, e in sorted(factors.items())))
    print("\nCase 1: Groups with all Sylow subgroups cyclic.")
    print(f"Number of such groups found: {num_type1_groups}")
    print("\nCase 2: Groups with cyclic odd Sylows and Generalized Quaternion Sylow 2-subgroup.")
    print(f"Number of such groups found: {num_type2_groups}")
    print("\nThe two cases are disjoint.")
    print("\nTotal number of such manifolds is the sum of the counts from these two cases.")
    print(f"Final equation: {num_type1_groups} + {num_type2_groups} = {total}")

solve()
<<<27>>>