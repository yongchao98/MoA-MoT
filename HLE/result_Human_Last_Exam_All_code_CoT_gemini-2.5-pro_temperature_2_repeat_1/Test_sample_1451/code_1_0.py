def solve():
    """
    Calculates T(n) = 2n - 1 for a given list of n values.
    """
    ns = [2, 3, 1234, 6712]
    results = []

    print("The formula for the minimum number of trials is T(n) = 2n - 1.")
    print("-" * 30)

    for n in ns:
        t_n = 2 * n - 1
        results.append(str(t_n))
        print(f"For n={n}: T({n}) = 2 * {n} - 1 = {t_n}")

    print("-" * 30)
    final_answer = ",".join(results)
    print(f"The values of T(2), T(3), T(1234), and T(6712) are:")
    print(final_answer)

solve()