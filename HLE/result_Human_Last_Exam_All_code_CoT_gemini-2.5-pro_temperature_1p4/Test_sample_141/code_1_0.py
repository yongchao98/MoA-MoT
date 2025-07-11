import math

def solve_a21():
    """
    Calculates the value of a_21 based on the derived recurrence relations.
    a_n is the number of ways to partition an n-element set into non-empty subsets,
    with a weak order on each subset.
    """
    
    N = 21

    # Step 1: Pre-compute binomial coefficients C(n, k)
    C = [[0] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        C[n][0] = 1
        for k in range(1, n + 1):
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k]

    # Step 2: Compute ordered Bell numbers f_k up to N
    # f_n = sum_{k=0 to n-1} C(n,k) * f_k
    f = [0] * (N + 2)
    f[0] = 1
    for n in range(1, N + 2):
        s = 0
        for k in range(n):
            s += C[n][k] * f[k]
        f[n] = s

    # Step 3: Compute a_n up to N
    # a_{n+1} = sum_{k=0 to n} C(n,k) * a_k * f_{n-k+1}
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(N):
        s = 0
        for k in range(n + 1):
            s += C[n][k] * a[k] * f[n - k + 1]
        a[n + 1] = s
        
    # Step 4: Compute and print the final equation for a_21
    # The recurrence for a_21 is for n=20:
    # a_21 = sum_{k=0 to 20} C(20,k) * a_k * f_{21-k}
    
    terms = []
    for k in range(N + 1):
        term = C[N-1][k] * a[k] * f[N - k]
        terms.append(term)
        
    print(f"The value of a_21 is found by the summation a_21 = sum_{k=0 to 20} C(20,k) * a_k * f_{21-k}.")
    print("The final equation is:")
    
    # Building the string for the equation
    equation_str = " + ".join([str(t) for t in terms])
    print(f"a_21 = {equation_str}")
    
    # Calculate the final sum
    final_value = sum(terms)
    print(f"\nThe exact numerical value of a_21 is:")
    print(f"a_21 = {final_value}")

solve_a21()