import math

def combinations(n, k):
    """Helper function to calculate combinations (n choose k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_part1():
    """
    Calculates the maximum number of exams for n=14.
    """
    print("### Part 1: Finding the maximum number of exams for n = 14 ###\n")
    
    n = 14
    k = 4
    
    print(f"Given parameters: n = {n} (total questions), k = {k} (questions per exam).\n")
    print("The core constraint is that any two exams share at most one question.")
    print("This implies that any pair of questions {q_a, q_b} can appear in at most one exam.\n")
    
    # Step 1: Find the maximum number of times a single question can appear (r_max)
    print("Step 1: Find the maximum number of exams any single question can be a part of.")
    print("Let 'r' be the number of exams a question appears in. For each such exam, the other k-1 questions must be distinct from the other questions in the other r-1 exams.")
    print(f"This gives the inequality: r * (k - 1) <= n - 1")
    k_minus_1 = k - 1
    n_minus_1 = n - 1
    print(f"Substituting values: r * {k_minus_1} <= {n_minus_1}")
    r_max = math.floor(n_minus_1 / k_minus_1)
    print(f"So, r <= {n_minus_1 / k_minus_1:.2f}, which means r_max = {r_max}.\n")

    # Step 2: Find the maximum number of exams (m)
    print("Step 2: Use r_max to find an upper bound for the number of exams 'm'.")
    print("The total number of question slots across all exams is m * k.")
    print("This must be less than or equal to the total possible slots, n * r_max.")
    print(f"The inequality is: m * k <= n * r_max")
    m_bound = math.floor((n * r_max) / k)
    print(f"Substituting values: m * {k} <= {n} * {r_max}")
    print(f"m * {k} <= {n * r_max}")
    print(f"So, m <= {n * r_max / k:.2f}, which means m <= {m_bound}.\n")

    # Step 3: Discuss the existence of the configuration for m = 14 and m = 13.
    print("Step 3: Check if a configuration with m = 14 is possible.")
    print(f"If m = 14, then the equality 4 * 14 = 14 * 4 must hold, which means every question must appear in exactly r = 4 exams.")
    print("A system of 14 exams on 14 questions, each of size 4, where every question appears 4 times and any two exams intersect at most once is a known combinatorial problem.")
    print("It has been proven that such a configuration does not exist.\n")
    
    print("Since m=14 is not possible, the maximum must be m <= 13.")
    print("A configuration for m=13 exams on n=13 questions exists (the projective plane of order 3, PG(2,3)).")
    print("We can use this configuration with 13 of our 14 questions, which satisfies all conditions.")
    
    answer_1 = 13
    print(f"\nConclusion for Part 1: The maximum number of different exams is {answer_1}.")
    return answer_1

def solve_part2():
    """
    Calculates the minimum n required for 10 exams.
    """
    print("\n### Part 2: Finding the minimum n to prepare 10 exams ###\n")
    m = 10
    k = 4
    
    print(f"Given parameters: m = {m} (number of exams), k = {k} (questions per exam).\n")
    
    # Step 1: Lower bound for n from counting pairs
    print("Step 1: Find a lower bound for 'n' by counting pairs of questions.")
    pairs_per_exam = combinations(k, 2)
    total_pairs_needed = m * pairs_per_exam
    print(f"Each exam has C({k}, 2) = {pairs_per_exam} pairs of questions.")
    print(f"For {m} exams, we need {m} * {pairs_per_exam} = {total_pairs_needed} distinct pairs of questions.")
    print("The total number of pairs available from 'n' questions is C(n, 2) = n*(n-1)/2.")
    print(f"So, we must have n*(n-1)/2 >= {total_pairs_needed}, which means n*(n-1) >= {2 * total_pairs_needed}.")
    
    n_min_pairs = 0
    while n_min_pairs * (n_min_pairs - 1) < 2 * total_pairs_needed:
        n_min_pairs += 1
    print(f"Testing values for n: 11*10 = 110 < 120. 12*11 = 132 >= 120.")
    print(f"The smallest integer 'n' that satisfies this is {n_min_pairs}.\n")

    # Step 2: Check feasibility for n >= n_min_pairs
    print("Step 2: Check if n = 12 is feasible using the single-question bound.")
    
    # Check n = 12
    n_test = 12
    total_incidences = m * k
    print(f"Total question slots to be filled = m * k = {m} * {k} = {total_incidences}.")
    print(f"This sum must equal the sum of appearances of each question (r_j): sum(r_j) = {total_incidences}.")
    
    print(f"\nChecking n = {n_test}:")
    n_minus_1 = n_test - 1
    k_minus_1 = k - 1
    r_max_n12 = math.floor(n_minus_1 / k_minus_1)
    print(f"The max appearances for any question is r_max <= ({n_test}-1)/({k}-1) = {n_minus_1/k_minus_1:.2f}, so r_max = {r_max_n12}.")
    max_total_incidences_n12 = n_test * r_max_n12
    print(f"The maximum possible sum of appearances is n * r_max = {n_test} * {r_max_n12} = {max_total_incidences_n12}.")
    print(f"This leads to a contradiction: the required sum is {total_incidences}, but the maximum possible is {max_total_incidences_n12}. So n=12 is impossible.\n")

    # Check n = 13
    n_test = 13
    print(f"Checking n = {n_test}:")
    n_minus_1 = n_test - 1
    r_max_n13 = math.floor(n_minus_1 / k_minus_1)
    print(f"The max appearances for any question is r_max <= ({n_test}-1)/({k}-1) = {n_minus_1/k_minus_1:.2f}, so r_max = {r_max_n13}.")
    max_total_incidences_n13 = n_test * r_max_n13
    print(f"The maximum possible sum of appearances is n * r_max = {n_test} * {r_max_n13} = {max_total_incidences_n13}.")
    print(f"This is consistent, since required sum {total_incidences} <= {max_total_incidences_n13}.")
    
    print("\nA configuration for n=13 and m=10 is known to exist. We can take 10 of the 13 exams from the projective plane PG(2,3).")
    
    answer_2 = 13
    print(f"\nConclusion for Part 2: The minimum value of n is {answer_2}.")
    return answer_2

if __name__ == '__main__':
    ans1 = solve_part1()
    ans2 = solve_part2()
    # The final answer format is not required here, but would be <<<[13, 13]>>> if needed
    