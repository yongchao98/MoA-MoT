import math

def solve_exam_problem():
    """
    Solves the two-part math competition exam problem using combinatorial principles.
    """
    # General parameters
    k = 4  # Questions per exam

    # --- Part 1: If n = 14, how many different exams can be created? ---
    print("--- Part 1: Find the maximum number of exams (m) for n=14 questions ---")
    n1 = 14
    print(f"Given: Number of questions n = {n1}, Questions per exam k = {k}\n")

    print("Step 1: Derive upper bounds for m using combinatorial inequalities.")
    
    # Bound 1: Based on counting pairs of questions
    print("\nBound 1: Each pair of questions can appear in at most one exam.")
    pairs_per_exam = math.comb(k, 2)
    total_pairs = math.comb(n1, 2)
    print(f"Number of pairs in one exam = C(k, 2) = C({k}, 2) = {pairs_per_exam}")
    print(f"Total possible pairs from n questions = C(n, 2) = C({n1}, 2) = {total_pairs}")
    print(f"The inequality is m * {pairs_per_exam} <= {total_pairs}")
    m_bound1 = total_pairs / pairs_per_exam
    print(f"m <= {total_pairs} / {pairs_per_exam} = {m_bound1:.2f}")
    print(f"Therefore, m <= {math.floor(m_bound1)}")

    # Bound 2: Based on question occurrences (Johnson Bound)
    print("\nBound 2: Based on how many exams one question can be in.")
    # For a question q, let r be the number of exams containing q.
    # The other k-1 questions in these r exams must be distinct.
    # r * (k-1) <= n-1
    max_r = (n1 - 1) // (k - 1)
    print(f"Let r be the number of exams any single question appears in.")
    print(f"r * (k-1) <= n-1  => r * {k-1} <= {n1-1} => r <= {(n1-1)/(k-1):.2f}")
    print(f"So, the maximum occurrences for any question is r_max = {max_r}.")
    # m * k = sum(r_j) <= n * r_max
    print(f"The inequality is m * k <= n * r_max => m * {k} <= {n1} * {max_r}")
    m_bound2 = (n1 * max_r) / k
    print(f"m <= ({n1} * {max_r}) / {k} = {m_bound2:.2f}")
    print(f"Therefore, m <= {math.floor(m_bound2)}")
    
    final_m_bound = math.floor(m_bound2)
    print(f"\nThe tighter bound is m <= {final_m_bound}.")

    print("\nStep 2: Check if m = 14 is possible.")
    print("If m = 14, every question must appear in exactly r = 4 exams for the bound to hold.")
    print("This would form a (v, k, r) = (14, 4, 4) symmetric design.")
    print("For such a design, every pair of questions must appear in lambda exams, where lambda = r(k-1)/(v-1).")
    r, v = 4, 14
    lambda_val = r * (k - 1) / (v - 1)
    print(f"lambda = {r}({k}-1) / ({v}-1) = {r * (k-1)} / {v-1} = {lambda_val:.4f}")
    print("Since lambda is not an integer, a design with m = 14 is not possible.")

    print("\nStep 3: Check if m = 13 is possible.")
    print("A known combinatorial object, the projective plane of order 3, is a set of 13 points and 13 lines (exams).")
    print("Each line contains 4 points, and any two lines intersect at exactly one point.")
    print("This provides a valid construction for m = 13 exams on n = 13 questions.")
    print("We can use this construction for n = 14 by simply not using the 14th question.")
    print("Thus, m = 13 is achievable.")
    
    ans1 = 13
    print("\nFinal equation for Part 1:")
    print(f"Maximum number of exams for n = {n1} is {ans1}")

    # --- Part 2: What is the minimum value of n needed to prepare 10 exams? ---
    print("\n\n--- Part 2: Find the minimum number of questions (n) for m=10 exams ---")
    m2 = 10
    print(f"Given: Number of exams m = {m2}, Questions per exam k = {k}\n")

    print("Step 1: Use the bounds to find a lower limit for n.")
    
    # Bound 1: Based on pairs
    print("\nBound 1: n(n-1) >= m * k(k-1)")
    min_prod1 = m2 * math.comb(k, 2)
    print(f"n(n-1) >= {m2} * {math.comb(k, 2)} = {min_prod1}")
    n_bound1 = 0
    n_val = k
    while True:
        if n_val * (n_val - 1) >= min_prod1:
            n_bound1 = n_val
            break
        n_val += 1
    print(f"Testing n = {n_bound1 - 1}: {n_bound1 - 1}*{n_bound1 - 2} = {(n_bound1-1)*(n_bound1-2)} < {min_prod1}")
    print(f"Testing n = {n_bound1}: {n_bound1}*{n_bound1-1} = {n_bound1*(n_bound1-1)} >= {min_prod1}")
    print(f"From this, n must be at least {n_bound1}.")
    
    # Bound 2: Based on occurrences
    print("\nBound 2: n * floor((n-1)/(k-1)) >= m * k")
    min_prod2 = m2 * k
    print(f"n * floor((n-1)/{k-1}) >= {m2} * {k} = {min_prod2}")
    n_bound2 = 0
    n_val = n_bound1 # Start checking from the previous bound
    while True:
        prod = n_val * ((n_val - 1) // (k-1))
        if prod >= min_prod2:
            n_bound2 = n_val
            break
        n_val += 1
    
    prev_n = n_bound2 - 1
    prev_prod = prev_n * ((prev_n - 1)//(k-1))
    print(f"Testing n = {prev_n}: {prev_n} * floor(({prev_n-1})/{k-1}) = {prev_prod} < {min_prod2}")
    curr_prod = n_bound2 * ((n_bound2-1)//(k-1))
    print(f"Testing n = {n_bound2}: {n_bound2} * floor(({n_bound2-1})/{k-1}) = {curr_prod} >= {min_prod2}")
    print(f"From this, n must be at least {n_bound2}.")

    ans2 = n_bound2
    print(f"\nThe tightest lower bound is n >= {ans2}.")

    print("\nStep 2: Check if n = 13 is possible for m = 10.")
    print("As established in Part 1, it's possible to create 13 exams with 13 questions.")
    print("Therefore, we can simply choose any 10 of these 13 exams.")
    print("Thus, m = 10 exams with n = 13 questions is achievable.")

    print("\nFinal equation for Part 2:")
    print(f"Minimum number of questions for m = {m2} is {ans2}")

solve_exam_problem()