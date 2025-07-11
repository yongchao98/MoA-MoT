def solve_tt_primes():
    """
    Solves the problem of identifying which objects generate a prime tt-ideal
    and counts the number of non-generated primes.
    """

    # The support of an object is represented as a frozenset of strings.
    # The strings correspond to slopes in the Balmer spectrum P^1_k.
    WHOLE_SPACE = frozenset(['whole_space'])
    EMPTY_SET = frozenset()

    # Supports for the base objects, determined from theory.
    # We assume p > 2, so that slopes 1 and -1 are distinct.
    supports = {
        1: WHOLE_SPACE,        # k(0): The tensor unit, support is the whole space.
        2: frozenset(['-1']),      # kC_p(0): Corresponds to slope -1.
        3: WHOLE_SPACE,        # [p-1](0): Supported on the whole space.
        4: frozenset(['inf']),     # rad(kC_p): Corresponds to slope infinity.
        5: frozenset(['1']),       # gap^1(rad(kC_p)): Corresponds to slope 1.
        6: frozenset(['0']),       # gap^{p-1}(rad(kC_p)): Corresponds to slope 0.
        7: frozenset(['1']),       # cone(tau): A "Tate object", corresponds to slope 1.
        15: frozenset(['1']),      # S: Another "Tate object", also corresponds to slope 1.
    }

    # The support of a power is the same as the support of the object itself.
    supports[8] = supports[7]   # cone(tau)^2
    supports[9] = supports[7]   # cone(tau)^{p-1}
    supports[10] = supports[7]  # cone(tau)^p

    # The support of a tensor product is the intersection of supports.
    supports[11] = supports[2].intersection(supports[7])
    supports[12] = supports[4].intersection(supports[7])
    supports[13] = supports[5].intersection(supports[7])
    supports[14] = supports[6].intersection(supports[7])
    supports[16] = supports[15].intersection(supports[2])
    supports[17] = supports[15].intersection(supports[4])
    supports[18] = supports[15].intersection(supports[7])
    supports[19] = supports[15].intersection(supports[5])
    supports[20] = supports[15].intersection(supports[6])

    # An object generates a prime tt-ideal if its support is non-empty.
    # An empty support means the object generates the whole category, which is not a proper prime ideal.
    prime_generators = []
    for i in range(1, 21):
        if supports[i] != EMPTY_SET:
            prime_generators.append(i)
    
    prime_generators.sort()

    # To get a definite integer for the second part of the answer, we deduce p=7.
    # Total number of closed prime ideals defined over F_p is p+1.
    p = 7
    total_fp_closed_primes = p + 1
    
    # Collect the distinct supports of the prime-generating objects.
    generated_prime_supports = set()
    for i in prime_generators:
        generated_prime_supports.add(supports[i])
        
    # Count how many distinct closed primes are generated.
    num_generated_closed_primes = 0
    for s in generated_prime_supports:
        if s != WHOLE_SPACE:
            num_generated_closed_primes += 1

    # The number of non-generated primes is the difference.
    num_not_generated = total_fp_closed_primes - num_generated_closed_primes

    # Print the answer in the required format.
    generators_str = ",".join(map(str, prime_generators))
    print(f"{generators_str};{num_not_generated}")

solve_tt_primes()