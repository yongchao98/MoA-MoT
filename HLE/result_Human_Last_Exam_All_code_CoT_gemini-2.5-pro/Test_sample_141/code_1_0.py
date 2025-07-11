import math

def solve():
    """
    Computes the value of a_21, where a_n is the number of ways to partition
    an n-element set and then place a weak order on each subset.
    """
    N = 21

    # Step 1: Compute Fubini numbers F_k up to F_N.
    # F_k is the number of weak orderings on a k-element set.
    # Recurrence: F_k = sum_{i=0 to k-1} C(k,i) * F_i, with F_0 = 1.
    # We need F_k up to k=21 for the final calculation.
    F = [0] * (N + 2)
    F[0] = 1
    for k in range(1, N + 2):
        s = 0
        for i in range(k):
            s += math.comb(k, i) * F[i]
        F[k] = s

    # Step 2: Compute a_n up to a_N.
    # a_n is the number of partitions of an n-set into lists.
    # Recurrence: a_{n+1} = sum_{k=0 to n} C(n,k) * a_k * F_{n-k+1}, with a_0 = 1.
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(N):  # Loop from n=0 to 20, which computes a[1] to a[21]
        s = 0
        for k in range(n + 1):
            s += math.comb(n, k) * a[k] * F[n - k + 1]
        a[n + 1] = s
    
    # Step 3: Print the final calculation and the result.
    # The calculation for a_21 is based on the recurrence for n=20:
    # a_21 = sum_{k=0 to 20} C(20, k) * a_k * F_{21-k}
    
    print("The value of a_21 is determined by the recurrence relation:")
    print("a_21 = sum_{k=0 to 20} C(20, k) * a_k * F_{21-k}")
    print("\nThe equation with the computed integer values for each term is:")
    
    terms = []
    # Loop from k=0 to 20 to calculate the 21 terms of the sum
    for k in range(N):
        # The term in the sum for a_21 (where n=20):
        # term_k = C(20, k) * a_k * F_{20-k+1} = C(20, k) * a_k * F_{21-k}
        term = math.comb(20, k) * a[k] * F[21 - k]
        terms.append(str(term))
        
    # Format the sum with line breaks for readability
    print("a_21 = " + " + \n       ".join(terms))
    
    final_result = a[N]
    print(f"\nThe final numerical value of a_21 is:")
    print(final_result)

solve()