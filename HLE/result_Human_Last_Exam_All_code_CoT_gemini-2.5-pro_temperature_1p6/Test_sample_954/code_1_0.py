def solve():
    """
    This function implements the plan to find the three-digit code.
    Based on the analysis, the maximum number of steps f(n) for a forest
    with maximum degree D is determined by a caterpillar graph construction.
    This leads to f(n) being on the order of D^2 when n is large enough,
    and O(n) otherwise.

    Let's break it down for each case.
    """

    # Case 1: Maximum degree D <= sqrt(log n).
    # For this case, D is very small compared to n. We can construct a caterpillar
    # graph that requires Theta(D^2) steps.
    # f_1(n) = Theta(D^2) = Theta((sqrt(log n))^2) = Theta(log n).
    # Looking at the categories provided:
    # 1) O(1)
    # 2) O(log^{0.1}log n)
    # 3) O(log^{0.9} log n)
    # 4) O(log log n)
    # 5) 2^{O(sqrt(log log n))}
    # 6) O(log^{0.9} n)
    # 7) o(log n)
    # 8) Theta(log n)
    # 9) omega(log n)
    # The complexity Theta(log n) matches category 8.
    f1_category = 8

    # Case 2: Maximum degree D <= log n.
    # Similar to case 1, D is small enough compared to n that we can achieve
    # the Theta(D^2) bound.
    # f_2(n) = Theta(D^2) = Theta((log n)^2).
    # We need to categorize Theta((log n)^2). This function grows faster
    # than log n. Therefore, it is omega(log n).
    # This matches category 9.
    f2_category = 9

    # Case 3: Any forest.
    # The maximum degree D can be large, up to n-1.
    # In the caterpillar construction, the number of steps is O(k*D) where k = min(D, n/D).
    # If D is large (e.g., D > sqrt(n)), then k is limited by n/D.
    # The number of steps is O((n/D) * D) = O(n).
    # The complexity f_3(n) = O(n).
    # O(n) grows faster than log n, so it is omega(log n).
    # This also matches category 9.
    f3_category = 9

    # Combine the digits to form the final three-digit number.
    final_answer = f"{f1_category}{f2_category}{f3_category}"
    print(final_answer)

solve()