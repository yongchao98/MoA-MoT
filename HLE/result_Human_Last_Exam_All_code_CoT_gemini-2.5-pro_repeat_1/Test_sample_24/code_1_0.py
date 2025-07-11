import math

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    INDEX = 7

    # n2[k] will store the number of elements sigma in S_k such that sigma^2 = id.
    # n5[k] will store the number of elements tau in S_k such that tau^5 = id.
    n2 = [0] * (INDEX + 1)
    n5 = [0] * (INDEX + 1)
    
    # N[k] will store the total number of homomorphisms from G to S_k.
    # T[k] will store the number of transitive homomorphisms from G to S_k.
    N = [0] * (INDEX + 1)
    T = [0] * (INDEX + 1)

    # Base cases
    n2[0] = 1
    n2[1] = 1
    for i in range(5):
      n5[i] = 1
    
    N[0] = 1

    print("Step 1: Calculate N_k (total homomorphisms to S_k) and T_k (transitive homomorphisms to S_k) for k = 1 to 7.")
    print("-" * 80)
    print(f"{'k':>2} | {'|s^2=1|':>8} | {'|t^5=1|':>8} | {'N_k':>10} | {'T_k':>10}")
    print("-" * 80)

    for k in range(1, INDEX + 1):
        # Calculate n2[k] using the recurrence: n2[k] = n2[k-1] + (k-1)*n2[k-2]
        if k >= 2:
            n2[k] = n2[k - 1] + (k - 1) * n2[k - 2]
        
        # Calculate n5[k] using the recurrence: n5[k] = n5[k-1] + (k-1)...(k-4)*n5[k-5]
        if k >= 5:
            p = 1
            for i in range(1, 5):
                p *= (k - i)
            n5[k] = n5[k - 1] + p * n5[k - 5]
        
        # N_k = n2[k] * n5[k]
        N[k] = n2[k] * n5[k]

        # Calculate T_k using the recurrence T_k = N_k - sum(C(k-1, i-1)*T_i*N_{k-i})
        intransitive_sum = 0
        for i in range(1, k):
            intransitive_sum += math.comb(k - 1, i - 1) * T[i] * N[k - i]
        
        T[k] = N[k] - intransitive_sum
        
        print(f"{k:>2} | {n2[k]:>8} | {n5[k]:>8} | {N[k]:>10} | {T[k]:>10}")

    print("-" * 80)
    
    # The number of subgroups of index n is T_n / (n-1)!
    t_n = T[INDEX]
    n_minus_1_factorial = math.factorial(INDEX - 1)
    num_subgroups = t_n // n_minus_1_factorial

    print("\nStep 2: Calculate the number of subgroups of index 7.")
    print(f"The number of transitive homomorphisms T_7 is {t_n}.")
    print(f"The number of subgroups of index 7 is given by T_7 / (7-1)!")
    print(f"Final Calculation: {t_n} / {n_minus_1_factorial} = {num_subgroups}")
    
    # Final answer in the required format
    print(f"\n<<<{num_subgroups}>>>")

solve()