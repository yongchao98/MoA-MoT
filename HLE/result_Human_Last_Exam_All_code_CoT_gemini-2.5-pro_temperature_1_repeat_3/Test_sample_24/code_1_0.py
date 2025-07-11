import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_i2(k):
    """Calculates the number of elements of order dividing 2 in S_k (involutions)."""
    # Recurrence relation: i2(k) = i2(k-1) + (k-1)*i2(k-2)
    # i2(0) = 1, i2(1) = 1
    if k == 0: return 1
    if k == 1: return 1
    a, b = 1, 1
    for i in range(2, k + 1):
        a, b = b, b + (i - 1) * a
    return b

def calculate_i5(k):
    """Calculates the number of elements of order dividing 5 in S_k."""
    # Such elements consist of only 1-cycles and 5-cycles.
    # Formula: i5(k) = sum_{j=0}^{floor(k/5)} k! / (j! * (k-5j)! * 5^j)
    total = 0
    for j in range(k // 5 + 1):
        term = math.factorial(k) // (math.factorial(j) * math.factorial(k - 5 * j) * (5 ** j))
        total += term
    return total

def solve_subgroup_problem():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7
    i2 = [0] * (N + 1)
    i5 = [0] * (N + 1)
    h = [0] * (N + 1)
    t = [0] * (N + 1)
    
    # h[0] is needed for the recurrence, h_0 represents a homomorphism to S_0 (trivial group)
    # i2[0] = 1, i5[0]=1 => h[0]=1
    h[0] = 1

    print("Step 1: Calculate i_2(k), i_5(k), and h_k = i_2(k) * i_5(k) for k=1 to 7.")
    for k in range(1, N + 1):
        i2[k] = calculate_i2(k)
        i5[k] = calculate_i5(k)
        h[k] = i2[k] * i5[k]
        print(f"k={k}: i_2({k})={i2[k]}, i_5({k})={i5[k]} => h_{k}={h[k]}")

    print("\nStep 2: Calculate t_k (number of transitive homomorphisms) for k=1 to 7.")
    for n in range(1, N + 1):
        sum_term = 0
        sum_str_parts = []
        for k in range(1, n):
            term = combinations(n - 1, k - 1) * t[k] * h[n - k]
            sum_term += term
            if term > 0:
                sum_str_parts.append(f"C({n-1},{k-1})*t_{k}*h_{n-k}")
        
        t[n] = h[n] - sum_term
        
        if n == 1:
            print(f"t_1 = h_1 = {h[1]}")
        else:
            print(f"\nt_{n} = h_{n} - ( {' + '.join(sum_str_parts)} )")
            sum_val_str = []
            current_sum = 0
            for k in range(1, n):
                C_nk = combinations(n-1, k-1)
                term = C_nk * t[k] * h[n-k]
                if term > 0:
                    sum_val_str.append(f"{C_nk}*{t[k]}*{h[n-k]}")
                current_sum += term
            print(f"t_{n} = {h[n]} - ( {' + '.join(sum_val_str)} )")
            print(f"t_{n} = {h[n]} - {current_sum}")
            print(f"t_{n} = {t[n]}")

    print("\nStep 3: Calculate the number of subgroups of index 7.")
    num_subgroups = t[N] // math.factorial(N - 1)
    print(f"Number of subgroups = t_7 / (7-1)!")
    print(f"Number of subgroups = {t[N]} / {math.factorial(N - 1)}")
    print(f"Number of subgroups = {num_subgroups}")
    
    return num_subgroups

final_answer = solve_subgroup_problem()
print(f"\n<<<{final_answer}>>>")
