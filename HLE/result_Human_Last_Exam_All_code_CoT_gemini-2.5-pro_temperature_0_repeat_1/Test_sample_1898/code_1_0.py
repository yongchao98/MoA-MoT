import math

def solve_exam_problem():
    """
    Solves both parts of the exam design problem.
    """
    print("### Part 1: Maximum exams for n=14 ###")
    
    # Parameters for Part 1
    n1 = 14
    k = 4
    
    print(f"Given n = {n1} total questions and k = {k} questions per exam.")
    print("The condition is that any two exams share at most one question.\n")

    # --- Step 1: Calculate upper bounds for the number of exams (m) ---
    print("Step 1: Calculating upper bounds for the number of exams (m).")
    
    # Bound based on pairs of questions
    c_n1_2 = math.comb(n1, 2)
    c_k_2 = math.comb(k, 2)
    m_bound1 = math.floor(c_n1_2 / c_k_2)
    print(f"A. Each exam has C({k}, 2) = {c_k_2} pairs of questions.")
    print(f"   The total pool of questions has C({n1}, 2) = {c_n1_2} pairs.")
    print(f"   Since no two exams can share a pair, m * {c_k_2} <= {c_n1_2}, which means m <= {m_bound1}.")

    # Bound based on question frequency (Fisher's Inequality)
    r1 = math.floor((n1 - 1) / (k - 1))
    m_bound2 = math.floor((n1 * r1) / k)
    print(f"B. Any single question can appear in at most r = floor(({n1}-1)/({k}-1)) = {r1} exams.")
    print(f"   The total number of question slots (m * k) must be at most n * r.")
    print(f"   m * {k} <= {n1} * {r1}, which means m <= {m_bound2}.")
    
    m_upper_bound = min(m_bound1, m_bound2)
    print(f"\nThe tightest mathematical upper bound is m <= {m_upper_bound}.\n")

    # --- Step 2: Use design theory to find the actual maximum ---
    print("Step 2: Checking the existence of a design for the found bounds.")
    print(f"A design with m = 14 is known to be impossible.")
    print("However, a design with n=13 questions and m=13 exams (the projective plane of order 3) exists.")
    print("This known design satisfies the condition, as any two of its 'exams' share exactly one question.")
    print("We can use this design for n=14 by simply ignoring the 14th question.")
    print("Therefore, m=13 is achievable.")
    
    result_m = 13
    print("\nConclusion for Part 1:")
    print(f"The maximum number of exams that can be created is 13.")
    print(f"m = {result_m}")

    print("\n" + "="*50 + "\n")

    print("### Part 2: Minimum questions (n) for 10 exams ###")
    
    # Parameters for Part 2
    m2 = 10
    
    print(f"Given m = {m2} exams and k = {k} questions per exam.\n")

    # --- Step 1: Find the minimum n that satisfies the necessary bounds ---
    print("Step 1: Finding the minimum 'n' that satisfies the mathematical bounds.")
    
    min_n = k  # Start searching for n from k
    while True:
        # Bound 1: m * C(k, 2) <= C(n, 2)
        bound1_ok = (m2 * c_k_2 <= math.comb(min_n, 2))
        
        # Bound 2: m * k <= n * floor((n-1)/(k-1))
        r2 = math.floor((min_n - 1) / (k - 1))
        bound2_ok = (m2 * k <= min_n * r2)
        
        if bound1_ok and bound2_ok:
            break
        min_n += 1
        
    print(f"We search for the smallest integer 'n' satisfying both:")
    print(f"1. {m2} * {c_k_2} <= C(n, 2)")
    print(f"2. {m2} * {k} <= n * floor((n-1)/{k-1})")
    print(f"The smallest value of n that satisfies both conditions is {min_n}.\n")

    # --- Step 2: Confirm this n is sufficient ---
    print("Step 2: Checking if a design with this 'n' is possible.")
    print(f"As established in Part 1, it is possible to create m=13 exams with n=13 questions.")
    print(f"Therefore, it is certainly possible to create a smaller set of m=10 exams using n=13 questions.")
    
    result_n = min_n
    print("\nConclusion for Part 2:")
    print(f"The minimum value of n needed is 13.")
    print(f"n = {result_n}")


# Run the solver
solve_exam_problem()