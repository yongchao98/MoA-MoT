import math
from collections import defaultdict

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n and their powers."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def euler_totient(n):
    """Calculates Euler's totient function phi(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return result

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def mobius_mu(n):
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def solve():
    """
    Finds the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Let d = {d} and g = {g}.")
    print("We want to find the number of primitive Dirichlet characters of conductor d and order g.")
    print("-" * 30)

    # Step 1: Factorize d
    d_factors = get_prime_factorization(d)
    prime_factors_d = sorted(d_factors.keys())
    print(f"Step 1: The prime factorization of d = {d} is {' * '.join(map(str, prime_factors_d))}.")
    print("-" * 30)

    # Step 2: Explain method and calculate Psi(k) for each k|g
    print("Step 2: We calculate Psi(k), the number of primitive characters with conductor d")
    print("and order dividing k, for each divisor k of g.")
    
    g_divisors = get_divisors(g)
    print(f"The divisors of g = {g} are {g_divisors}.")
    
    Psi_values = {}
    for k in g_divisors:
        print(f"\n--- Calculating Psi({k}) ---")
        
        prod_M_p_k = 1
        M_p_k_expressions = []
        
        for p in prime_factors_d:
            p_minus_1 = p - 1
            k_divisors = get_divisors(k)
            
            sum_phi_j = 0
            j_list = []
            phi_j_values = []

            for j in k_divisors:
                if j > 1 and p_minus_1 % j == 0:
                    j_list.append(j)
                    phi_val = euler_totient(j)
                    phi_j_values.append(phi_val)
                    sum_phi_j += phi_val
            
            M_p_k = sum_phi_j
            print(f"For p = {p}, p-1 = {p_minus_1}. The orders j must divide k={k}, be greater than 1, and divide {p_minus_1}.")
            if not j_list:
                print(f"  No such j > 1 exists. So M_{p}({k}) = 0.")
            else:
                print(f"  The possible orders j are {j_list}.")
                phi_expr = ' + '.join([f"phi({j})" for j in j_list])
                val_expr = ' + '.join(map(str, phi_j_values))
                print(f"  M_{p}({k}) = {phi_expr} = {val_expr} = {M_p_k}.")
            
            prod_M_p_k *= M_p_k
            M_p_k_expressions.append(str(M_p_k))
            
        Psi_values[k] = prod_M_p_k
        print(f"Psi({k}) = {' * '.join(M_p_k_expressions)} = {prod_M_p_k}.")

    print("-" * 30)
    # Step 3: Calculate the final answer using Mobius inversion
    print("Step 3: We use the Mobius inversion formula to find the number of characters of order exactly g.")
    print(f"psi(g) = sum_{{k|g}} mu(g/k) * Psi(k)")
    
    final_answer = 0
    terms = []
    full_terms = []
    
    for k in g_divisors:
        mu_val = mobius_mu(g // k)
        psi_val = Psi_values[k]
        term = mu_val * psi_val
        final_answer += term
        terms.append(f"({mu_val}) * ({psi_val})")
        full_terms.append(f"mu({g//k})*Psi({k})")

    print(f"psi({g}) = {' + '.join(full_terms)}")
    print(f"psi({g}) = {' + '.join(terms)}")
    print(f"psi({g}) = {final_answer}")
    print("-" * 30)
    
    print(f"The number of primitive Dirichlet characters of conductor {d} and order {g} is {final_answer}.")

solve()