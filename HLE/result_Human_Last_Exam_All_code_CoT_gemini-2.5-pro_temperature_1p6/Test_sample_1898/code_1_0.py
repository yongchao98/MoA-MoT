import math

def solve_exam_problem():
    """
    Solves the two-part exam preparation problem using a step-by-step analytical approach.
    The code will print the reasoning and calculations for each part.
    """

    print("This program solves a combinatorial problem about creating exams.")
    print("Each exam has 4 questions, and any two exams can share at most one question.\n")

    # --- Part 1: If n = 14, how many different exams can be created? ---
    print("--- Part 1: Maximum number of exams for n = 14 questions ---\n")
    n1 = 14
    k = 4

    print(f"Let n be the number of questions ({n1}), and k be the number of questions per exam ({k}).")
    print("Let m be the maximum number of exams we want to find.\n")

    print("Step 1: Find an upper bound for m.")
    print("Let's consider a single question q. Let it appear in r_q exams.")
    print(f"For each of these r_q exams, the other k-1 = {k-1} questions must be distinct from the questions in the other exams containing q.")
    print(f"These {k-1}-question sets are disjoint and drawn from the remaining n-1 = {n1-1} questions.")
    print(f"Therefore, the total number of distinct questions needed is r_q * (k-1) <= n-1.")
    
    r_max = math.floor((n1 - 1) / (k - 1))
    print(f"This gives an upper bound on r_q: r_q <= floor((n-1)/(k-1)) = floor(({n1}-1)/({k}-1)) = {r_max}.\n")

    print("The total number of question 'slots' across all m exams is m * k.")
    print("This is also the sum of occurrences for each question: Sum(r_q) <= n * r_max.")
    print(f"So, m * {k} <= {n1} * {r_max}, which implies m <= floor(({n1} * {r_max}) / {k}).")
    
    m_upper_bound = math.floor((n1 * r_max) / k)
    print(f"Equation: m <= floor(({n1} * {r_max}) / {k}) = {m_upper_bound}\n")
    print(f"So, the maximum number of exams is at most {m_upper_bound}.\n")
    
    print("Step 2: Establish a lower bound for m.")
    print("For n=13, it's a known result in combinatorial design theory that a structure called a 'projective plane of order 3' exists.")
    print("This structure consists of 13 exams on 13 questions, where every pair of questions appears in exactly one exam.")
    print("This implies any two exams share exactly one question, satisfying our condition.")
    print(f"Since 13 exams can be made with 13 questions, at least 13 can be made with {n1} questions (by not using the 14th question).")
    print(f"So, we know that m >= 13. This means m is either 13 or 14.\n")

    print("Step 3: Prove by contradiction that m cannot be 14.")
    print(f"Assume for contradiction that m = {m_upper_bound} is possible.")
    print(f"From Step 1, this would mean every question must appear in exactly r = {r_max} exams.")
    
    total_exam_pairs = math.comb(m_upper_bound, 2)
    total_intersections = n1 * math.comb(r_max, 2)
    disjoint_pairs = total_exam_pairs - total_intersections

    print(f"Let's count the total number of intersections between all pairs of exams: Sum(|Ei intersect Ej|).")
    print(f"From the question perspective, each of the {n1} questions appears in {r_max} exams, contributing to C({r_max},2) = {math.comb(r_max, 2)} intersections.")
    print(f"Total Intersections = {n1} * {math.comb(r_max, 2)} = {total_intersections}.")
    print(f"Since |Ei intersect Ej| is either 0 or 1, this sum means {total_intersections} pairs of exams share one question.")
    print(f"The total number of pairs of exams is C(m,2) = C({m_upper_bound},2) = {total_exam_pairs}.")
    print(f"This implies there are {total_exam_pairs} - {total_intersections} = {disjoint_pairs} pairs of exams that are disjoint (share no questions).\n")

    print("If m=14, there must be 7 disjoint pairs of exams. Let (E1, E2) be one such pair.")
    print(f"Let E1 = {{q1, q2, q3, q4}}. It is disjoint from E2. E1 must intersect the other 12 exams (E3 to E14).")
    print(f"The 4 questions in E1 must each appear in r-1 = {r_max-1} other exams. These {r_max-1} exams for q1 must be different from the {r_max-1} exams for q2, etc.")
    print(f"This partitions the 12 exams (E3..E14) into 4 groups of 3, one group for each question in E1.")
    
    print("\nNow, consider exam E2. Each of the 12 exams (E3 to E14) must also intersect E2 in exactly one question.")
    print(f"Consider one group of 3 exams, all containing q1. Each of these 3 exams must intersect E2 at a distinct point.")
    print("(If two shared an intersection point p from E2, they would both contain {{q1, p}}, which is not allowed).")
    print(f"So, each of the 4 groups of exams requires {r_max-1} distinct questions from E2 to serve as intersection points.")
    required_qs_from_e2 = 4 * (k-1)
    print(f"In total, this requires 4 groups * {k-1} questions/group = {required_qs_from_e2} distinct questions from E2.")
    print(f"But E2 only contains k = {k} questions. It's impossible to provide {required_qs_from_e2} distinct questions from a set of {k}.")
    print("This is a contradiction. The assumption that m=14 is false.\n")

    q1_answer = 13
    print(f"Conclusion for Part 1: Since m >= 13 and m < 14, the maximum number of exams is {q1_answer}.\n")


    # --- Part 2: What is the minimum value of n needed to prepare 10 exams? ---
    print("--- Part 2: Minimum number of questions n for m = 10 exams ---\n")
    m2 = 10

    print("Step 1: Use the inequality from Part 1 to find a lower bound for n.")
    print(f"We need to find the smallest n such that m <= floor((n/k) * floor((n-1)/(k-1))).")
    print(f"We check for n until {m2} <= floor((n/{k}) * floor((n-1)/{k-1})).\n")

    print("Testing values for n:")
    n2 = k # Start checking from n=k
    while True:
        n2 += 1
        inner_floor = math.floor((n2 - 1) / (k - 1))
        m_bound = math.floor(n2 / k * inner_floor)
        
        # We need each part of the equation printed
        n_div_k_str = f"{n2}/{k}"
        nm1_div_km1_str = f"({n2}-1)/({k}-1)"
        print(f"For n = {n2}: m <= floor(({n_div_k_str}) * floor({nm1_div_km1_str})) = floor({n2/k:.2f} * {inner_floor}) = {m_bound}")
        
        if m_bound >= m2:
            q2_answer = n2
            print(f"At n={n2}, the bound is {m_bound}, which is >= {m2}. This is the first n that could work.\n")
            break
        else:
            print(f"m_bound ({m_bound}) is less than {m2}. So n={n2} is not enough.")
            
    print(f"Step 2: Check if n = {q2_answer} is sufficient.")
    print(f"The calculation shows we need at least n = {q2_answer} questions.")
    print(f"As established in Part 1, for n={q2_answer}, it's possible to create up to {m_bound} exams.")
    print(f"Since we need only {m2} exams, and {m2} <= {m_bound}, this is achievable.\n")

    print(f"Conclusion for Part 2: The minimum value of n needed is {q2_answer}.\n")

if __name__ == '__main__':
    solve_exam_problem()