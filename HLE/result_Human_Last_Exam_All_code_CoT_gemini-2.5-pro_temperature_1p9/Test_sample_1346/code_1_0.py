def solve_and_print():
    """
    Calculates the required values based on the logic described above.
    The code calculates a(n) iteratively for the required small values of n
    and shows the calculation steps as requested.
    """
    results = []

    # For p = 50051, we need to calculate a(5).
    p1 = 50051
    n1 = 5
    print(f"# For p = {p1}, we determined the value is a({n1}).")
    a_prev = 1
    a_curr = 3
    print("a(0) = 1")
    print("a(1) = 3")
    for i in range(2, n1 + 1):
        a_next = 4 * a_curr - a_prev
        print(f"a({i}) = 4 * {a_curr} - {a_prev} = {a_next}")
        a_prev = a_curr
        a_curr = a_next
    result1 = a_curr
    results.append(result1)
    print(f"# The result for p={p1} is {result1}.\n")

    # For p = 50069, we need to calculate a(3).
    p2 = 50069
    n2 = 3
    print(f"# For p = {p2}, we determined the value is a({n2}).")
    a_prev = 1
    a_curr = 3
    print("a(0) = 1")
    print("a(1) = 3")
    for i in range(2, n2 + 1):
        a_next = 4 * a_curr - a_prev
        print(f"a({i}) = 4 * {a_curr} - {a_prev} = {a_next}")
        a_prev = a_curr
        a_curr = a_next
    result2 = a_curr
    results.append(result2)
    print(f"# The result for p={p2} is {result2}.\n")

    # Print the final comma-separated result.
    final_answer = ",".join(map(str, results))
    print(final_answer)

solve_and_print()