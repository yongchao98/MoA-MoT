import math

def prime_factorize(n):
    """Returns the list of prime factors of a number."""
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.append(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    return factors

def get_divisors(n):
    """Returns the set of divisors of a number."""
    divs = {1}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    divs.add(n)
    return sorted(list(divs))

def phi(n):
    """Calculates Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def solve():
    """
    Finds the number of primitive Dirichlet characters of a given conductor d and order g.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.\n")

    # Step 1: Prime factorize the conductor d
    p_factors = prime_factorize(d)
    m = len(p_factors)
    print(f"Step 1: The prime factorization of d = {d} is {p_factors}.")
    print(f"The number of prime factors is m = {m}.\n")

    # Step 2: For this conductor (product of distinct primes), the number of primitive characters
    # of order g is given by the inclusion-exclusion formula:
    # N(d, g) = N_g(d) - N_{g/p1}(d) - N_{g/p2}(d) - ... + N_{g/(p1*p2)}(d) + ...
    # For g=6, primes are 2,3. So, N(d,6) = N_6(d) - N_3(d) - N_2(d) + N_1(d).
    print("Step 2: We use the inclusion-exclusion principle for order g = 6:")
    print("Result = N_6(d) - N_3(d) - N_2(d) + N_1(d)\n")
    
    # Step 3: Calculate N_k(d) for k in {1, 2, 3, 6}
    print("Step 3: Calculating each term N_k(d).")
    # N_k(d) is the number of primitive characters with conductor d and order dividing k.
    # It is calculated as (C_k)^m, where C_k is the number of non-trivial characters mod p
    # with order dividing k, and m is the number of prime factors of d.
    # This assumes C_k is the same for all prime factors, which we will verify.
    
    memo_C = {}

    for k in [6, 3, 2, 1]:
        divs_k = [div for div in get_divisors(k) if div > 1]
        
        num_choices_per_prime = 0
        all_primes_compatible = True
        
        # Check compatibility and calculate C_k
        first_prime = True
        for j in divs_k:
            is_compatible_for_all_p = True
            for p in p_factors:
                if (p - 1) % j != 0:
                    is_compatible_for_all_p = False
                    break
            if not is_compatible_for_all_p:
                all_primes_compatible = False
                print(f"For order j={j}, not all p-1 are divisible by it. This case needs special handling (not required here).")
                break
            num_choices_per_prime += phi(j)
        
        if not all_primes_compatible:
            # This problem has compatible primes, so this branch won't be hit.
            memo_C[k] = -1 # Sentinel for error
            continue

        Nk_d = num_choices_per_prime**m
        memo_C[k] = (num_choices_per_prime, Nk_d)
        
        phi_calcs = " + ".join([f"phi({j})" for j in divs_k]) if divs_k else "0"
        phi_vals = " + ".join([str(phi(j)) for j in divs_k]) if divs_k else "0"
        print(f"For k = {k}:")
        print(f"  Non-trivial orders j dividing k are {divs_k}.")
        print(f"  Number of choices per prime C({k}) = {phi_calcs} = {phi_vals} = {num_choices_per_prime}.")
        print(f"  N_{k}(d) = C({k})^m = {num_choices_per_prime}^{m} = {Nk_d}.")

    # Final calculation
    N6_d = memo_C[6][1]
    N3_d = memo_C[3][1]
    N2_d = memo_C[2][1]
    N1_d = memo_C[1][1]

    result = N6_d - N3_d - N2_d + N1_d
    
    print("\nStep 4: Final calculation.")
    print(f"The number of primitive Dirichlet characters is N_6(d) - N_3(d) - N_2(d) + N_1(d)")
    print(f"= {N6_d} - {N3_d} - {N2_d} + {N1_d}")
    print(f"= {result}")
    
    return result

final_answer = solve()
print(f"\nFinal Answer: {final_answer}")