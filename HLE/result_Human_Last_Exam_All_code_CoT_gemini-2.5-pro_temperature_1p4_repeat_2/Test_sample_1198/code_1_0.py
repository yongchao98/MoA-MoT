import math

def search_a_mod_2():
    """
    This function searches for a real number 'a' such that
    floor(a^n) is congruent to n mod 2 for all n > 0.
    It works by iteratively narrowing an interval [low, high) for 'a'.
    """
    # Start with a reasonably wide interval for a. Let's start with floor(a) = 1.
    low = 1.0
    high = 2.0
    
    # We will check for n from 1 up to N
    N = 20

    print("--- Searching for 'a' in the case modulo 2 ---")
    print(f"Initial interval for a: [{low:.10f}, {high:.10f})")

    possible = True
    for n in range(1, N + 1):
        # Determine the required parity for floor(a^n)
        target_parity = n % 2
        
        # Current range for a^n is [low^n, high^n)
        a_n_low = low**n
        a_n_high = high**n

        # We need to find an integer k in the range of a^n with the correct parity.
        # k = floor(a^n).
        
        # Start searching for a suitable integer k.
        # The smallest possible k is floor(a_n_low).
        k_start = math.floor(a_n_low)
        
        found_k_for_n = False
        # Iterate through candidate values for k. We only need to find one that works.
        # The range of k can be large, but the smallest candidates are most likely to work.
        for k_candidate in range(k_start, math.ceil(a_n_high) + 5):
            if k_candidate > 0 and k_candidate % 2 == target_parity:
                # This k_candidate has the correct parity.
                # It implies a new interval for 'a': [k^(1/n), (k+1)^(1/n))
                k_low = k_candidate**(1/n)
                k_high = (k_candidate + 1)**(1/n)
                
                # We need to intersect this new interval with our existing [low, high)
                new_low = max(low, k_low)
                new_high = min(high, k_high)
                
                # Use a small epsilon for floating point comparisons
                epsilon = 1e-15
                if new_low < new_high - epsilon:
                    # We found a valid choice for k. Update the interval for 'a'.
                    low = new_low
                    high = new_high
                    print(f"n={n:2}: Chose k={k_candidate}. New interval for a: [{low:.10f}, {high:.10f})")
                    found_k_for_n = True
                    break # Success for this n, move to the next.
        
        if not found_k_for_n:
            print(f"n={n:2}: Could not find a suitable k. The interval has become empty.")
            possible = False
            break

    print("\nConclusion:")
    if possible:
        print("The search was successful up to N={}. This provides strong evidence that".format(N))
        print("such a number 'a' exists for the modulo 2 case.")
    else:
        print("The search failed. This suggests no such 'a' exists or a different choice of k was needed.")

search_a_mod_2()
