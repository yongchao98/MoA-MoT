import math

def find_counterexample_mod_3():
    """
    Tests for the existence of a real number a > 0 such that
    floor(a^n) is congruent to n (mod 3) for all integers n > 0.
    
    A known mathematical result states that no such 'a' exists.
    A step in the proof shows that if 'a' existed, a^3 must be an integer
    multiple of 3. We test candidates of this form, a = (3k)^(1/3),
    to find a contradiction.
    """
    print("Investigating if there is an a > 0 such that floor(a^n) % 3 == n % 3 for all n > 0.")
    print("We test candidates of the form a = (3*k)^(1/3) for a contradiction.")
    print("-" * 30)

    # We will use k=1, which gives a = 3^(1/3)
    k = 1
    a = 3**(1/3)
    max_n = 10

    print(f"Testing candidate a = (3*k)^(1/3) for k={k}")
    print(f"a = 3^(1/3) ≈ {a:.5f}")
    
    for n in range(1, max_n + 1):
        # Calculate floor(a^n)
        val = math.floor(a**n)
        
        # Get the required and actual remainders modulo 3
        required_rem = n % 3
        actual_rem = val % 3
        
        print(f"n={n}: floor(a^{n}) = {val}, val % 3 = {actual_rem}. Required: {required_rem}", end="")
        
        if actual_rem == required_rem:
            print(" -> OK")
        else:
            print(" -> FAIL")
            print("-" * 30)
            print(f"A contradiction is found for a = 3^(1/3) at n = {n}.")
            return

find_counterexample_mod_3()