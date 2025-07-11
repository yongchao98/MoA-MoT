import math

def solve_exam_problem():
    """
    Solves both parts of the exam design problem and prints the step-by-step reasoning.
    """
    # --- Part 1: Find maximum m for n=14 ---
    print("Part 1: Finding the maximum number of exams (m) for n=14 questions.")
    print("----------------------------------------------------------------------")
    n1 = 14
    k = 4

    # Step 1: Bound from question pairs
    pairs_per_exam = math.comb(k, 2)
    total_pairs = math.comb(n1, 2)
    m_bound1 = total_pairs // pairs_per_exam
    
    print("Step 1: Calculate the bound based on pairs of questions.")
    print(f"Each exam has k={k} questions, containing C(k, 2) = C({k}, 2) = {pairs_per_exam} pairs.")
    print(f"The total pool of n={n1} questions has C(n, 2) = C({n1}, 2) = {total_pairs} available pairs.")
    print("The sets of pairs from each exam must be disjoint. This gives the inequality:")
    print(f"m * {pairs_per_exam} <= {total_pairs}")
    print(f"m <= {total_pairs / pairs_per_exam:.2f}")
    print(f"As m must be an integer, this implies m <= {m_bound1}.\n")

    # Step 2: Bound from question frequency
    r_max = (n1 - 1) // (k - 1)
    m_bound2 = (n1 * r_max) // k

    print("Step 2: Calculate a tighter bound based on individual question frequency.")
    print(f"Any single question can appear in at most r_max = floor((n-1)/(k-1)) = floor(({n1}-1)/({k}-1)) = {r_max} exams.")
    print("The total number of question instances across all exams (m * k) cannot exceed the maximum possible sum of frequencies (n * r_max).")
    print(f"m * {k} <= {n1} * {r_max}")
    print(f"m <= ({n1} * {r_max}) / {k} = {n1 * r_max / k}")
    print(f"This implies m <= {m_bound2}.\n")

    # Step 3: Conclusion for Part 1
    print("Step 3: Conclusion for Part 1.")
    print("The bounds show m <= 14. However, the existence of a configuration for m=14 is a known impossibility in combinatorial design theory.")
    print("The next highest value is 13. A configuration for m=13 exams with n=13 questions (a Steiner System S(2,4,13)) is known to exist.")
    print("We can use this known design and add a 14th question to the pool that is not used in any exam.")
    print("This creates a valid set of 13 exams for n=14 questions.")
    final_answer_part1 = 13
    print(f"\nFinal Answer for Part 1: The maximum number of different exams is {final_answer_part1}.")
    print("=" * 70 + "\n")

    # --- Part 2: Find minimum n for m=10 ---
    print("Part 2: Finding the minimum number of questions (n) for m=10 exams.")
    print("----------------------------------------------------------------------")
    m2 = 10

    print("Step 1: Find the minimum n that satisfies the necessary conditions (the two bounds).")
    print(f"We require n such that n(n-1) >= m*k*(k-1) = {m2}*{k}*({k-1}) = {m2 * k * (k - 1)}")
    print(f"and n * floor((n-1)/(k-1)) >= m*k = {m2}*{k} = {m2 * k}")
    
    # Search for the smallest n >= k that satisfies these conditions.
    min_n_candidate = k
    while True:
        cond1 = min_n_candidate * (min_n_candidate - 1) >= m2 * k * (k - 1)
        cond2 = min_n_candidate * ((min_n_candidate - 1) // (k - 1)) >= m2 * k
        
        print(f"\nTesting n = {min_n_candidate}:")
        val1 = min_n_candidate * (min_n_candidate - 1)
        print(f"  Condition 1: {val1} >= {m2 * k * (k - 1)}? {'Yes' if cond1 else 'No'}")
        val2 = min_n_candidate * ((min_n_candidate - 1) // (k - 1))
        print(f"  Condition 2: {val2} >= {m2 * k}? {'Yes' if cond2 else 'No'}")

        if cond1 and cond2:
            break
        min_n_candidate += 1
    
    final_answer_part2 = min_n_candidate
    print(f"\nBoth conditions are first met for n = {final_answer_part2}.\n")
    
    # Step 2: Conclusion for Part 2
    print("Step 2: Conclusion for Part 2.")
    print(f"The necessary conditions require that n must be at least {final_answer_part2}.")
    print("To confirm this is sufficient, we note that a set of 13 exams can be constructed with 13 questions.")
    print("Since it is possible to create 13 exams, it is certainly possible to prepare 10 of them.")
    print(f"\nFinal Answer for Part 2: The minimum number of questions needed is {final_answer_part2}.")


if __name__ == '__main__':
    solve_exam_problem()