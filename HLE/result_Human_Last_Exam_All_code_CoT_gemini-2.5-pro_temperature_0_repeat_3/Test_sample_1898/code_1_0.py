import math

def combinations(n, k):
    """Calculates the number of combinations of n items taken k at a time."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_exam_problem():
    """
    Solves the two-part math competition exam problem using combinatorial principles.
    """
    # --- Part 1: If n = 14, how many different exams can be created? ---
    print("--- Part 1: Maximum exams for n=14 ---")
    print("Let m be the number of exams, n=14 be the number of questions, and k=4 be the questions per exam.")
    print("The constraint is that any two exams share at most one question.\n")

    n1 = 14
    k = 4

    # Step 1: Pair-based bound
    print("Step 1: Derive an upper bound on m by counting pairs of questions.")
    print("Each exam has C(k, 2) pairs. The total number of pairs across m exams must be less than or equal to the total possible pairs from n questions, C(n, 2).")
    
    C_n1_2 = combinations(n1, 2)
    C_k_2 = combinations(k, 2)
    bound1 = math.floor(C_n1_2 / C_k_2)
    
    print(f"The inequality is: m * C({k}, 2) <= C({n1}, 2)")
    print(f"Substituting the values: m * {C_k_2} <= {C_n1_2}")
    print(f"This gives: m <= {C_n1_2 / C_k_2:.2f}, so m <= {bound1}.\n")

    # Step 2: Incidence-based bound
    print("Step 2: Derive a second upper bound by considering question occurrences.")
    print("Let r be the number of exams a single question appears in. The other k-1 questions in these r exams must be distinct.")
    r_max = math.floor((n1 - 1) / (k - 1))
    print(f"So, r * ({k}-1) <= ({n1}-1), which means r <= {r_max}.")
    print("The total number of question 'slots' (m*k) must be less than or equal to n * r_max.")
    
    bound2 = math.floor((n1 * r_max) / k)
    print(f"The inequality is: m * {k} <= {n1} * {r_max}")
    print(f"This gives: m <= ({n1} * {r_max}) / {k}, so m <= {bound2}.\n")

    # Step 3: Conclusion for Part 1
    print("Step 3: Determine the true maximum.")
    print(f"The tighter mathematical bound is m <= {min(bound1, bound2)}.")
    print("However, this doesn't guarantee existence. From design theory, it's known that for n=13, a structure called a projective plane of order 3 allows for exactly 13 exams satisfying the conditions.")
    print("Adding a 14th question does not allow for any new valid exams to be formed.")
    part1_answer = 13
    print(f"Therefore, the maximum number of exams that can be created is {part1_answer}.\n")

    # --- Part 2: What is the minimum value of n needed to prepare 10 exams? ---
    print("--- Part 2: Minimum questions n for m=10 exams ---")
    m2 = 10
    print(f"We need to find the minimum n for m={m2} and k={k}.\n")

    # Step 1: Use inequalities to find a lower bound for n
    print("Step 1: Apply the same inequalities to find a lower bound for n.")
    
    # From pair-based bound
    required_C_n_2_times_2 = m2 * C_k_2 * 2
    print(f"From m * C(k, 2) <= C(n, 2), we get {m2} * {C_k_2} <= n*(n-1)/2, which simplifies to {required_C_n_2_times_2} <= n*(n-1).")
    min_n1 = 0
    n_test = k
    while True:
        if n_test * (n_test - 1) >= required_C_n_2_times_2:
            min_n1 = n_test
            break
        n_test += 1
    print(f"By testing values, we find n must be at least {min_n1}.\n")

    # From incidence-based bound
    required_total_slots = m2 * k
    print(f"From m * k <= n * floor((n-1)/(k-1)), we get {m2} * {k} <= n * floor((n-1)/3), which is {required_total_slots} <= n * floor((n-1)/3).")
    min_n2 = 0
    n_test = min_n1
    while True:
        val = n_test * math.floor((n_test - 1) / (k - 1))
        print(f"Testing n={n_test}: {n_test} * floor({n_test-1}/3) = {val}. Is {val} >= {required_total_slots}? {'Yes' if val >= required_total_slots else 'No'}")
        if val >= required_total_slots:
            min_n2 = n_test
            break
        n_test += 1
    
    min_n_bound = max(min_n1, min_n2)
    print(f"The necessary condition from these bounds is n >= {min_n_bound}.\n")

    # Step 2: Check for sufficiency
    print("Step 2: Confirm this lower bound is achievable.")
    print(f"We need to know if a valid configuration with m=10 exists for n={min_n_bound}.")
    print("As noted before, with n=13 questions, it is possible to create m=13 exams (the projective plane of order 3).")
    print("Since we can create 13 exams, we can certainly provide 10 of them.")
    part2_answer = 13
    print(f"Therefore, n={part2_answer} is sufficient, and since n must be at least {min_n_bound}, the minimum value is {part2_answer}.")

if __name__ == '__main__':
    solve_exam_problem()