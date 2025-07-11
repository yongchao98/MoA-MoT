import math

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    MAX_N = 7

    # a2[n] = number of elements in S_n with order dividing 2 (involutions)
    # Recurrence: a2[n] = a2[n-1] + (n-1)*a2[n-2]
    a2 = [0] * (MAX_N + 1)
    a2[0] = 1
    a2[1] = 1
    for n in range(2, MAX_N + 1):
        a2[n] = a2[n-1] + (n-1) * a2[n-2]

    # a5[n] = number of elements in S_n with order dividing 5
    # Recurrence: a5[n] = a5[n-1] + (n-1)(n-2)(n-3)(n-4)*a5[n-5]
    a5 = [0] * (MAX_N + 1)
    for n in range(MAX_N + 1):
        if n < 5:
            a5[n] = 1
        else:
            term = math.perm(n - 1, 4)
            a5[n] = a5[n-1] + term * a5[n-5]

    # h[n] = |Hom(G, S_n)| = a2[n] * a5[n]
    h = [0] * (MAX_N + 1)
    for n in range(MAX_N + 1):
        h[n] = a2[n] * a5[n]

    # N[n] = number of subgroups of index n
    N = [0] * (MAX_N + 1)
    
    # Calculate N[n] recursively
    for n in range(1, MAX_N + 1):
        if n == 1:
            N[n] = h[n]
            continue
        
        sum_term = 0
        for k in range(1, n):
            # This is N_k * (n-1)!/(n-k)! * h_{n-k}
            term = N[k] * math.perm(n - 1, k - 1) * math.factorial(k-1) * h[n-k]
            sum_term += term
        
        N[n] = (h[n] - sum_term) // math.factorial(n - 1)

    # Print the detailed calculation for N_7
    n = 7
    print(f"The number of subgroups of index {n} in G = C_2 * C_5 is N_{n}.")
    print("This is calculated using the formula:")
    print(f"N_{n} * ({n-1})! = h_{n} - sum_{{k=1}}^{{{n-1}}} [N_k * ({n-1})!/({n}-k)! * h_{{{n}-k}}]\n")

    print("First, we list the required values:")
    for i in range(1, n + 1):
        print(f"N_{i} = {N[i]}")
    print("-" * 20)
    for i in range(n + 1):
        print(f"h_{i} = |Hom(G, S_{i})| = {h[i]}")
    print("-" * 20)

    print(f"Plugging the values into the formula for N_{n}:")
    
    sum_terms_str = []
    sum_terms_val = []
    for k in range(1, n):
        perm_val = math.perm(n - 1, k)
        # The formula uses (n-1)!/(n-k)! which is perm(n-1, k)
        # N_k * (n-1)!/(n-k)! * h_{n-k}
        term_val = N[k] * math.perm(n-1, k) * h[n-k]
        sum_terms_val.append(term_val)
        sum_terms_str.append(f"({N[k]} * {perm_val} * {h[n-k]})")

    total_sum = sum(sum_terms_val)
    lhs_val = N[n] * math.factorial(n-1)

    print(f"N_{n} * {n-1}! = {h[n]} - [ " + " + ".join(sum_terms_str) + " ]")
    print(f"{N[n]} * {math.factorial(n-1)} = {h[n]} - [ " + " + ".join(map(str, sum_terms_val)) + " ]")
    print(f"{lhs_val} = {h[n]} - {total_sum}")
    print(f"{lhs_val} = {h[n] - total_sum}")
    print(f"\nTherefore, N_{n} = {lhs_val} / {math.factorial(n-1)} = {N[n]}")
    
    print("\nThe number of subgroups of index 7 is:")
    print(N[n])


solve()