import math

def solve_exam_problem():
    """
    Solves the two-part exam creation problem by providing a step-by-step
    analytical explanation and calculating the final answers.
    """
    k = 4  # questions per exam

    # --- Part 1: n = 14, find max m ---
    n1 = 14
    print("--- Part 1: Maximum exams for n=14 ---")
    print(f"Given n = {n1} questions, find the maximum number of exams (m).")
    print(f"Each exam has k = {k} questions.")
    print("Condition: Any two exams have at most one common question.\n")

    print("Step 1: Find the maximum number of times a single question can appear (r).")
    print("Consider a single question. Let it appear in 'r' exams.")
    print("Each of these 'r' exams contains this question and 3 others.")
    print("These sets of 3 other questions must be disjoint, otherwise two exams would share more than one question.")
    print(f"The total number of distinct questions in these sets is r * (k-1), which must be at most n-1.")
    r_max_val = math.floor((n1 - 1) / (k - 1))
    print(f"Equation: r * ({k}-1) <= {n1}-1  =>  r * 3 <= 13  =>  r <= floor(13/3) = {r_max_val}")
    print(f"So, a single question can appear in at most {r_max_val} exams.\n")

    print("Step 2: Use the limit on 'r' to find an upper bound for 'm'.")
    print("The total number of question 'slots' across all 'm' exams is m * k.")
    print("This must equal the sum of appearances for all 'n' questions.")
    print(f"Equation: m * k <= n * r_max  =>  m * {k} <= {n1} * {r_max_val}")
    m_bound = math.floor((n1 * r_max_val) / k)
    print(f"m <= ({n1} * {r_max_val}) / {k}  =>  m <= {m_bound}\n")
    
    print("Step 3: Analyze the problem using known combinatorial structures.")
    print("The bound m <= 14 is not guaranteed to be achievable.")
    print("It's a known result that for n=13, a system of 13 exams (called a Projective Plane) can be created.")
    print("Let's see if we can add a new exam by using the 14th question.")
    print("A new exam would need 3 questions from the original 13 and the new 14th question.")
    print("Let this new exam be E_new = {q1, q2, q3, 14}.")
    print("However, in the n=13 system, the pair of questions {q1, q2} already exists in some exam, say E_j.")
    print("This would cause E_new and E_j to share at least two questions, violating the main condition.")
    print("This shows we cannot simply add a new exam. The maximum is limited by the structure on 13 points.")
    
    answer_part1 = 13
    print(f"\nConclusion for Part 1: The maximum number of different exams that can be created is {answer_part1}.\n")

    # --- Part 2: m = 10, find min n ---
    m2 = 10
    print("--- Part 2: Minimum questions for m=10 ---")
    print(f"Given m = {m2} exams, find the minimum number of questions (n).\n")

    print("Step 1: Use the inequality from Part 1.")
    print(f"Inequality: m * k <= n * floor((n-1)/(k-1))")
    print(f"{m2} * {k} <= n * floor((n-1)/3)")
    print(f"40 <= n * floor((n-1)/3)\n")

    print("Step 2: Test integer values of 'n' to find the minimum that satisfies the inequality.")
    
    # Check n=12
    n_test_12 = 12
    val_12 = n_test_12 * math.floor((n_test_12 - 1) / 3)
    print(f"Checking for n = {n_test_12}:")
    print(f"Equation: {n_test_12} * floor(({n_test_12}-1)/3) = {n_test_12} * floor(11/3) = {n_test_12} * 3 = {val_12}")
    print(f"Is 40 <= {val_12}? {40 <= val_12}. So, n=12 is not sufficient.\n")

    # Check n=13
    n_test_13 = 13
    val_13 = n_test_13 * math.floor((n_test_13 - 1) / 3)
    print(f"Checking for n = {n_test_13}:")
    print(f"Equation: {n_test_13} * floor(({n_test_13}-1)/3) = {n_test_13} * floor(12/3) = {n_test_13} * 4 = {val_13}")
    print(f"Is 40 <= {val_13}? {40 <= val_13}. So, n=13 is a possible minimum.\n")

    print("Step 3: Confirm that a solution for n=13 exists.")
    print("As established in Part 1, it's possible to create 13 exams with 13 questions.")
    print("Therefore, to prepare 10 exams, we can simply use 10 of those 13, and n=13 is sufficient.")
    
    answer_part2 = 13
    print(f"\nConclusion for Part 2: The minimum value of n needed is {answer_part2}.")

solve_exam_problem()