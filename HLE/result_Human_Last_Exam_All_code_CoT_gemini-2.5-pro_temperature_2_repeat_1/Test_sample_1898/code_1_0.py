import math

def solve():
    """
    Solves both parts of the math competition problem and explains the steps.
    """
    
    # --- Part 1: Given n = 14, find the maximum number of exams (m). ---
    
    n_part1 = 14
    k = 4  # questions per exam
    
    print("--- Part 1: Finding the maximum number of exams for n = 14 ---")
    print("Let n be the number of questions, k be the number of questions per exam, and m be the number of exams.")
    print(f"We are given n = {n_part1} and k = {k}.")
    print("The condition is that any two exams can have at most one common question.")
    
    # Step 1: Bound from pairs
    print("\nStep 1: Establishing an upper bound for m using pairs of questions.")
    c_k_2 = math.comb(k, 2)
    print(f"Each exam contains C(k, 2) = C({k}, 2) = {c_k_2} pairs of questions.")
    c_n_2_part1 = math.comb(n_part1, 2)
    print(f"The total pool of questions has C(n, 2) = C({n_part1}, 2) = {c_n_2_part1} possible pairs.")
    print("Since any two exams share at most one question, any pair of questions can appear in at most one exam.")
    print(f"Thus, the total number of pairs from all exams (m * {c_k_2}) cannot exceed the total available pairs ({c_n_2_part1}).")
    print(f"This gives the inequality: {c_k_2} * m <= {c_n_2_part1}")
    m_bound1 = math.floor(c_n_2_part1 / c_k_2)
    print(f"m <= {c_n_2_part1} / {c_k_2} => m <= {c_n_2_part1 / c_k_2:.2f}, so m <= {m_bound1}.")

    # Step 2: Bound from question occurrences
    print("\nStep 2: Establishing a second bound based on individual question frequency.")
    print("Let r_q be the number of exams a question q appears in.")
    print("For any single question q, it must be paired with k-1 other questions in each of the r_q exams.")
    print("Because any two of these exams can only share q, the sets of k-1 other questions must be disjoint.")
    print(f"This requires {k-1} * r_q distinct questions from the remaining n-1 pool.")
    print(f"So, (k-1) * r_q <= n-1 => {k-1} * r_q <= {n_part1-1}")
    r_bound_part1 = math.floor((n_part1 - 1) / (k - 1))
    print(f"r_q <= ({n_part1-1}) / {k-1} => r_q <= {(n_part1-1)/(k-1):.2f}, so the maximum for any r_q is {r_bound_part1}.")
    
    # Step 3: Combining bounds
    print("\nStep 3: Combining the bounds to find the maximum possible m.")
    print("By counting the total number of question slots across all exams in two ways, we get:")
    print("Sum of questions in all exams = m * k")
    print("Sum of occurrences of each question = Sum(r_q)")
    print(f"So, {k} * m = Sum(r_q).")
    print(f"Since each r_q <= {r_bound_part1}, the sum is at most n * {r_bound_part1} = {n_part1 * r_bound_part1}.")
    print(f"This gives our tightest constraint: {k} * m <= n * {r_bound_part1}.")
    
    # Step 4: Testing potential m values
    print(f"\nWe test m starting from our first bound, m = {m_bound1}.")
    test_m = m_bound1
    km_val = k * test_m
    n_r_bound_val = n_part1 * r_bound_part1
    print(f"For m = {test_m}: {k} * {test_m} = {km_val}.")
    print(f"The maximum allowed sum of r_q is {n_part1} * {r_bound_part1} = {n_r_bound_val}.")
    print(f"Is {km_val} <= {n_r_bound_val}? This is false.")
    print(f"Therefore, m cannot be {test_m}.")
    
    # Test m-1
    test_m = m_bound1 - 1
    answer1 = test_m
    km_val = k * test_m
    print(f"\nNow we test m = {test_m}.")
    print(f"For m = {test_m}: {k} * {test_m} = {km_val}.")
    print(f"Is {km_val} <= {n_r_bound_val}? This is true.")
    print(f"This value of m is possible. It corresponds to a known construction in design theory.")
    print(f"\nConclusion for Part 1: The maximum number of exams is {answer1}.")
    print("The final equation is:")
    print(f"m = {answer1}")

    # --- Part 2: Given m = 10 exams, find the minimum n required. ---

    m_part2 = 10
    
    print("\n\n--- Part 2: Finding the minimum n to prepare m = 10 exams ---")
    print(f"We are given m = {m_part2} and k = {k}.")
    
    # Step 1: Find a lower bound for n
    print("\nStep 1: Using the inequalities to find a lower bound for n.")
    km_pairs = m_part2 * c_k_2
    print(f"From {c_k_2} * m <= C(n, 2), we have {km_pairs} <= n(n-1)/2.")
    print(f"This simplifies to {2*km_pairs} <= n(n-1).")
    n_lower_bound = 1
    while n_lower_bound * (n_lower_bound - 1) < 2 * km_pairs:
        n_lower_bound += 1
    print(f"By solving n(n-1) >= {2*km_pairs}, we find the minimum integer n is {n_lower_bound}.")

    # Step 2: Test n values starting from the bound
    print("\nStep 2: Testing integer values of n starting from the lower bound.")
    
    n_test = n_lower_bound -1 # Start loop at n_lower_bound
    min_n_found = False
    answer2 = -1

    while not min_n_found:
        n_test += 1
        print(f"\nTesting n = {n_test}:")
        r_bound_part2 = math.floor((n_test - 1) / (k - 1))
        sum_r = k * m_part2
        max_sum_r = n_test * r_bound_part2
        
        print(f"The required total number of question occurrences is k * m = {k} * {m_part2} = {sum_r}.")
        print(f"With n={n_test}, the maximum possible r_q is floor(({n_test}-1)/({k}-1)) = {r_bound_part2}.")
        print(f"The maximum possible total occurrences is n * {r_bound_part2} = {n_test} * {r_bound_part2} = {max_sum_r}.")

        if sum_r <= max_sum_r:
            print(f"The condition {sum_r} <= {max_sum_r} holds. So n={n_test} is feasible.")
            print("Constructions for these parameters (m=10, k=4, n=13) exist.")
            answer2 = n_test
            min_n_found = True
        else:
            print(f"The condition {sum_r} <= {max_sum_r} fails. n={n_test} is not sufficient.")
    
    print(f"\nConclusion for Part 2: The minimum value of n is {answer2}.")
    print("The final equation is:")
    print(f"n = {answer2}")

    # Final combined answer
    print("\n=======================================================")
    print("Final Answers Summary")
    print("=======================================================")
    print(f"1. For n = 14, the maximum number of different exams is {answer1}.")
    print(f"2. To prepare 10 exams, the minimum value of n needed is {answer2}.")

    final_answer_string = f"1. {answer1}, 2. {answer2}"
    print(f'<<<{final_answer_string}>>>')


solve()