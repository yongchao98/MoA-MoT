import math

def solve():
    """
    Solves the user's request based on the step-by-step plan.
    """
    
    # Question 1: What is the expected value of X_19?
    n_19 = 19
    # For odd n, E[X_n] = (n^2 - 1) / 2
    x_19 = (n_19**2 - 1) / 2
    
    # Question 2: What is the expected value of X_20?
    # For even n, the distance between gifts is always odd.
    # The game ends when the distance is 2 or n-2. For even n, these are even.
    # The target distance is unreachable, so the expected time is infinite.
    x_20 = float('inf')
    
    # Question 3: In general, for odd n > 1, what is the expected value of X_n?
    # The formula is (n^2 - 1) / 2.
    # The output string for the formula.
    formula_n = "(n^2 - 1) / 2"

    # Question 4: For odd n > 30, what is the expected number of times that between the two gifts there are exactly 10 and n-12 friends.
    # This means distance is 11.
    # By symmetry and properties of random walks on cycles, the expected number of visits is n+1.
    visits_formula = "n + 1"
    
    # Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n > 1 value?
    # Yes, because the target distance n-2 is reachable for odd n.
    finite_time_answer = "yes"
    
    # To be clear in the output, let's represent infinity with a string.
    # and combine all answers in a comma-separated string.
    # The requested output for question 3 requires showing the equation. So let's write it out.
    n_squared_minus_1_div_2_str = "(n*n - 1) / 2" # for textual representation
    # The requested output for question 4 also needs to be in an equation form.
    n_plus_1_str = "n + 1"
    
    # The final output needs each number in the equation printed. Let's do that for the first case n=19.
    final_eq_19 = f"({n_19}*{n_19} - 1) / 2 = {int(x_19)}"
    
    
    print(f"{int(x_19)}, inf, {formula_n}, {visits_formula}, {finite_time_answer}")
    
solve()