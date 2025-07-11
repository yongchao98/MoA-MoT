import math

# Memoization caches
h_memo = {}
u_memo = {}
t_memo = {}
FACTORIALS = [math.factorial(i) for i in range(11)]

def count_elements_of_order_dividing_p(p, n):
    """
    Counts the number of permutations in S_n whose order divides p.
    This works when p is prime. A permutation has order dividing p if and only if
    all its cycles have length 1 or p.
    """
    if (p, n) in h_memo:
        return h_memo[(p, n)]

    count = 0
    # A permutation has cycle structure (1)^k (p)^l where k + p*l = n
    for l in range(n // p + 1):
        k = n - p * l
        # Number of permutations with this cycle structure:
        # n! / (k! * l! * 1^k * p^l)
        term = FACTORIALS[n] // (FACTORIALS[k] * FACTORIALS[l] * (p ** l))
        count += term
    
    h_memo[(p, n)] = count
    return count

def u(n):
    """
    Calculates u_n = |Hom(G, S_n)| for G = C_2 * C_5.
    u_n = |Hom(C_2, S_n)| * |Hom(C_5, S_n)|.
    """
    if n in u_memo:
        return u_memo[n]
    if n == 0:
        return 1
        
    h2_n = count_elements_of_order_dividing_p(2, n)
    h5_n = count_elements_of_order_dividing_p(5, n)
    result = h2_n * h5_n
    u_memo[n] = result
    return result

def t(n):
    """
    Calculates t_n, the number of transitive homomorphisms from G to S_n.
    Uses the recurrence: t_n = u_n - sum_{k=1 to n-1} C(n-1,k-1) * t_k * u_{n-k}
    """
    if n in t_memo:
        return t_memo[n]
    if n == 0:
        # No transitive actions on 0 elements by convention
        return 0
    if n == 1:
        # The only action on 1 element is transitive.
        # |Hom(G, S_1)| = 1, so t_1=1
        t_memo[1] = u(1)
        return u(1)

    sum_val = 0
    for k in range(1, n):
        term = math.comb(n - 1, k - 1) * t(k) * u(n - k)
        sum_val += term
    
    result = u(n) - sum_val
    t_memo[n] = result
    return result

def solve():
    """
    Calculates the number of subgroups of index 7 and prints the process.
    """
    N = 7
    print(f"Finding the number of subgroups of index {N} in G = C2 * C5.\n")
    print("Let a_n be the number of subgroups of index n.")
    print("Let u_n = |Hom(G, S_n)| and t_n = |transitive Hom(G, S_n)|.")
    print("We use the formula: a_n = t_n / (n-1)!\n")
    
    # Pre-compute all values up to N
    for i in range(1, N + 1):
        u(i)
        t(i)

    print("Calculated values for u_n and t_n:")
    print("n | u_n      | t_n")
    print("--|----------|---------")
    for i in range(1, N + 1):
        print(f"{i} | {u(i):<8} | {t(i)}")
    print("\n")

    # Final calculation for a_7
    t7 = t(N)
    fact6 = FACTORIALS[N - 1]

    print(f"The number of subgroups of index {N} is a_{N} = t_{N} / ({N}-1)! = t_{N} / {N-1}!\n")
    print(f"To find t_{N}, we use the formula:")
    print(f"t_{N} = u_{N} - sum_{{k=1}}^{{{N-1}}} (C({N-1}, k-1) * t_k * u_{{{N}-k}})\n")

    print(f"First, u_{N} = |Hom(C2, S_{N})| * |Hom(C5, S_{N})|")
    h2_7 = count_elements_of_order_dividing_p(2, 7)
    h5_7 = count_elements_of_order_dividing_p(5, 7)
    print(f"|Hom(C2, S7)| = {h2_7}")
    print(f"|Hom(C5, S7)| = {h5_7}")
    print(f"u_7 = {h2_7} * {h5_7} = {u(7)}\n")

    print(f"Now for the sum in the expression for t_7:")
    sum_expr = []
    for k in range(1, N):
        C = math.comb(N-1, k-1)
        tk = t(k)
        u_nk = u(N-k)
        sum_expr.append(f"C(6, {k-1})*t_{k}*u_{N-k} = {C}*{tk}*{u_nk}")
    print("Sum = " + " + ".join(sum_expr))

    sum_terms = []
    total_sum = 0
    for k in range(1, N):
        C = math.comb(N-1, k-1)
        tk = t(k)
        u_nk = u(N-k)
        term_val = C * tk * u_nk
        sum_terms.append(str(term_val))
        total_sum += term_val

    print("Sum = " + " + ".join(sum_terms))
    print(f"Sum = {total_sum}\n")
    
    print(f"So, t_7 = u_7 - Sum = {u(7)} - {total_sum} = {t7}\n")

    print(f"Finally, a_7 = t_7 / 6! = {t7} / {fact6}")
    
    result = t7 // fact6
    print(f"a_7 = {result}\n")
    print(f"The number of subgroups of index 7 is {result}.")
    print(f"<<<{result}>>>")

solve()
