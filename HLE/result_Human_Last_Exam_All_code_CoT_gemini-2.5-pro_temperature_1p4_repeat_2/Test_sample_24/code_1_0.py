import math

# Use memoization for efficiency
memo_i2 = {0: 1, 1: 1}
memo_a = {0: 1}
memo_t = {}

def factorial(n):
    return math.factorial(n)

def combinations(n, k):
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))

def i2(n):
    """Calculates the number of elements of order dividing 2 in S_n.
    Uses the recurrence i2(n) = i2(n-1) + (n-1)*i2(n-2).
    """
    if n not in memo_i2:
        memo_i2[n] = i2(n - 1) + (n - 1) * i2(n - 2)
    return memo_i2[n]

def i5(n):
    """Calculates the number of elements of order dividing 5 in S_n.
    For n < 10, these are just the identity and 5-cycles.
    """
    if n < 5:
        return 1
    # Number of 5-cycles is C(n, 5) * 4!
    num_5_cycles = combinations(n, 5) * factorial(4)
    return 1 + num_5_cycles

def a(n):
    """Calculates the total number of homomorphisms from G to S_n."""
    if n not in memo_a:
        memo_a[n] = i2(n) * i5(n)
    return memo_a[n]

def t(n):
    """Calculates the number of transitive homomorphisms from G to S_n."""
    if n not in memo_t:
        if n == 0:
            return 0
        # t_n = a_n - sum_{k=1 to n-1} C(n-1, k-1) * t_k * a_{n-k}
        sum_val = 0
        for k in range(1, n):
            term = combinations(n - 1, k - 1) * t(k) * a(n - k)
            sum_val += term
        
        result = a(n) - sum_val
        memo_t[n] = result
    return memo_t[n]

def solve():
    """Main function to perform the calculation and print the result."""
    n_target = 7
    
    # Pre-calculate all necessary values
    for k in range(1, n_target + 1):
        a(k)
        t(k)

    print(f"To find the number of subgroups of index {n_target}, we first find t_{n_target}, the number of transitive homomorphisms from G to S_{n_target}.")
    print(f"The formula is: t_n = a_n - sum_{{k=1 to n-1}} C(n-1, k-1) * t_k * a_{{n-k}}")
    print("\nFor n = 7, the equation is:")
    print("t_7 = a_7 - [ C(6,0)t_1*a_6 + C(6,1)t_2*a_5 + C(6,2)t_3*a_4 + C(6,3)t_4*a_3 + C(6,4)t_5*a_2 + C(6,5)t_6*a_1 ]")
    
    # Building the substituted equation string
    substitutions = []
    total_sum = 0
    for k in range(1, n_target):
        comb = combinations(n_target - 1, k - 1)
        tk = t(k)
        ank = a(n_target - k)
        term = comb * tk * ank
        total_sum += term
        substitutions.append(f"{comb}*{tk}*{ank}")
    
    print("\nSubstituting the calculated values:")
    print(f"t_7 = {a(n_target)} - [ {' + '.join(substitutions)} ]")

    # Building the summed term string
    sum_terms = []
    for k in range(1, n_target):
        sum_terms.append(str(combinations(n_target - 1, k - 1) * t(k) * a(n_target - k)))

    print("\nCalculating the terms in the sum:")
    print(f"t_7 = {a(n_target)} - [ {' + '.join(sum_terms)} ]")
    print(f"t_7 = {a(n_target)} - {total_sum}")
    val_t7 = t(n_target)
    print(f"t_7 = {val_t7}")

    print("\nFinally, the number of subgroups N_7 is t_7 divided by (7-1)!:")
    val_fact = factorial(n_target - 1)
    result = val_t7 // val_fact
    print(f"N_7 = t_7 / 6! = {val_t7} / {val_fact} = {result}")
    
    print(f"\nThus, the number of subgroups of index 7 in G = C_2 * C_5 is {result}.")
    print(f"<<<{result}>>>")

solve()
