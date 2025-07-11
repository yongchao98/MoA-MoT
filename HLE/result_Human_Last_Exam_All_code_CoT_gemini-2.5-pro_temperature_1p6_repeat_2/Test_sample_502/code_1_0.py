import math

def find_max_dimension():
    """
    Finds the largest possible dimension for R/I based on the theory of
    coinvariant rings and complex reflection groups.
    """
    group_order = 10000

    print("Step 1: Theoretical background")
    print(f"The dimension of the coinvariant ring R/I is bounded by the order of the acting group H, dim(R/I) <= |H|.")
    print("Equality holds if H is a complex reflection group.")
    print(f"H is a quotient of G, so |H| must divide |G| = {group_order}.")
    print("Goal: Find the largest divisor of 10000 that is the order of a complex reflection group in dimension n <= 10.\n")

    # Get all divisors of the group order in descending order
    divs = []
    for i in range(1, int(math.sqrt(group_order)) + 1):
        if group_order % i == 0:
            divs.append(i)
            if i*i != group_order:
                divs.append(group_order//i)
    divs.sort(reverse=True)

    print(f"Step 2: Searching for the maximum order k among divisors of {group_order}.")

    max_k = 0
    found_params = {}

    for k in divs:
        if max_k > 0: # Already found the largest possible k
            break
            
        print(f"\nChecking potential order k = {k}...")
        
        # We only need to check the G(m,p,n) family.
        # Max n is limited because n! grows very fast.
        for n in range(1, 11):
            try:
                n_fact = math.factorial(n)
            except ValueError: # n is too large
                n_fact = float('inf')

            if n_fact > k:
                # If n! > k, then k = m^n * n! / p is impossible for m, p >= 1
                # unless n_fact does not divide k properly, but it gets more complex.
                # A simpler argument: k*p = m^n*n!, so n! must divide k*p.
                # Since k is getting small, this loop can terminate early.
                if n >= 3: # we already know 3! has a factor of 3 not in 10000.
                    # print(f"  For n={n}, n!={n_fact} is not a factor of any divisor of 10000. Skipping.")
                    continue


            # Order equation: k = m^n * n! / p  => m^n / p = k / n!
            # Let Q = k / n!. Q must be rational, but let's check for integer Q first for simplicity
            if k % n_fact != 0:
                continue

            Q = k // n_fact
            
            # We need to find integers m, p such that m^n = p * Q and p|m.
            # Let's search for m. The condition m|pQ implies m^(n-1)|Q.
            # So m <= Q^(1/(n-1)) for n>1.
            limit_m = 0
            if n == 1:
                # m = pQ. Try small p.
                limit_m = k + 1 # Search for m up to k (e.g., p=1 => m=k)
            else: # n > 1
                try:
                    limit_m = int(Q**(1/(n-1))) + 2
                except (ValueError, OverflowError): # Q might be too large
                    limit_m = 10000 # Heuristic limit

            for m in range(1, limit_m):
                if (m**n) % Q == 0:
                    p = (m**n) // Q
                    if p > 0 and m % p == 0:
                        print(f"  SUCCESS: Found a reflection group for n={n} of order k={k}.")
                        print(f"  Group G({m},{p},{n}) is a valid reflection group.")
                        
                        max_k = k
                        found_params = {'k': k, 'n': n, 'm': m, 'p': p}
                        break # Found for this n, go to next k
            if max_k > 0:
                break # Found for this k, stop searching n
        if max_k > 0:
            break # Found highest k, stop searching k
    
    print("\nStep 3: Conclusion")
    if max_k > 0:
        k = found_params['k']
        n = found_params['n']
        m = found_params['m']
        p = found_params['p']
        
        print(f"The largest possible dimension is the largest order found, which is {max_k}.")
        print("This can be achieved, for example, with a group G of order 10000 that acts on C^10")
        print(f"via a representation whose image is the reflection group H = G(m={m}, p={p}, n={n}).")
        print("The order of this group H is calculated as:")
        n_fact = math.factorial(n)
        print(f"|G({m},{p},{n})| = ({m}^{n} * {n}!) / {p} = ({m**n} * {n_fact}) / {p} = {(m**n * n_fact) // p}")

        final_dim = max_k
    else:
        print("Could not find a suitable reflection group. The theoretical maximum may not be achievable.")
        final_dim = "Not found"

    print("\nFinal Answer:")
    print(f"The largest possible dimension for R/I is {final_dim}.")
    return final_dim

# Execute the function to get the result
result = find_max_dimension()
print(f'<<<{result}>>>')
