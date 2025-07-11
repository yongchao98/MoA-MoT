import math

def solve():
    """
    This function calculates and prints the answers to the user's questions.
    """
    # Question 1: What is the expected value of X_19?
    n_19 = 19
    # For odd n, the formula is (n^2 - 1) / 2
    ex_19 = (n_19**2 - 1) / 2
    
    # Question 2: What is the expected value of X_20?
    # For even n, the expected value is infinite.
    ex_20 = float('inf')

    # Question 3: In general, for odd n>1 what is the expected value of X_n?
    # The formula is (n^2 - 1) / 2. We represent this as a string.
    formula_odd_n = "(n**2 - 1) / 2"
    
    # Question 4: For odd n>30 what is the expected number of times that between the two gifts there are exactly 10 and n-12 friends.
    # This corresponds to a distance of 11 between the gifts. The number of friends is 10.
    # This is a classic result for certain random walks, where the expected number of visits to a state k levels away is k.
    # The number of friends is 10.
    exp_visits = 10
    
    # Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value?
    # Yes, because the expected time is finite for odd n.
    finite_time_prob_one = "Yes"

    # We need to output the formula for E_n for odd n literally in the final equation as per the prompt instructions.
    # I will print the formula text and its value for n=19 for clarity.
    
    # Outputting the results in a single comma-separated line.
    # We will print the computed values and the formula string.
    # The prompt "Remember in the final code you still need to output each number in the final equation!" is slightly ambiguous.
    # I will interpret it as printing the answers in order.
    # Let's re-read the format requirement. "Give these answers separated by a comma."
    
    answer_1 = int(ex_19)
    answer_2 = "inf" if ex_20 == float('inf') else int(ex_20)
    # The question for odd n is about the formula.
    # I will output the expression string. "(n^2-1)/2"
    answer_3 = "(n^2-1)/2" # Representing the formula as requested by the problem's nature.
    answer_4 = exp_visits
    answer_5 = finite_time_prob_one
    
    print(f"{answer_1},{answer_2},{answer_3},{answer_4},{answer_5}")

solve()