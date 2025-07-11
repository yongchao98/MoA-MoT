import math

def solve_and_print_answers():
    """
    This function calculates and prints the answers to the user's questions.
    """
    
    # Part 1: E(X_19)
    n_19 = 19
    ex_19_val = 2 * (n_19 - 1)
    
    # Part 2: E(X_20)
    ex_20_val = math.inf
    
    # Part 3: E(X_n) for odd n > 1
    ex_n_odd_formula = "2*(n-1)"
    
    # Part 4: Expected number of times for specific distance
    expected_visits = 2
    
    # Part 5: Finite time probability
    finite_time_prob_one = "Yes"
    
    # Print the answers as requested
    print(f"The expected value of X_19 is given by the equation: 2 * (19 - 1) = {ex_19_val}")
    print(f"The expected value of X_20 is: {ex_20_val}")
    print(f"The expected value of X_n for odd n > 1 is: {ex_n_odd_formula}")
    print(f"For odd n > 30, the expected number of times the distance between gifts is 11 is: {expected_visits}")
    print(f"It is true that the game ends in finite time with probability 1 for every odd n > 1: {finite_time_prob_one}")
    
    # This is the final answer string for the platform.
    # final_answer_string = f"{ex_19_val}, {ex_20_val}, {ex_n_odd_formula}, {expected_visits}, {finite_time_prob_one}"
    # print(f"\nComma-separated answers: {final_answer_string}")

solve_and_print_answers()