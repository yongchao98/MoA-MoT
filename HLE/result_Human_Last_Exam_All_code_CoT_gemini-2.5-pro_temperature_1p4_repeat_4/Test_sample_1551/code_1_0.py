import math

def solve_puzzle():
    """
    Solves the puzzle by calculating the requested values and printing them.
    """

    # --- Part 1: E[X_19] ---
    n_19 = 19
    # For an odd n, the expected time E[X_n] = (n^2 - 1) / 2
    ex_19 = (n_19**2 - 1) / 2

    # --- Part 2: E[X_20] ---
    n_20 = 20
    # For an even n, the parity of the distance between gifts is fixed and different from 0.
    # The game never ends.
    ex_20 = float('inf')

    # --- Part 3: E[X_n] for odd n > 1 ---
    # The general formula for odd n.
    general_formula_odd_n = "(n^2 - 1) / 2"

    # --- Part 4: Expected number of times for a specific configuration ---
    # The state of the system is the difference D_t. D_0 = n-1 is even for odd n.
    # The distance between gifts being 10 corresponds to D_t=11 or D_t=n-11.
    # D_t is always even, so it can never be 11. It can be n-11 (which is even).
    # The expected number of visits to an intermediate state for this type of random walk is 2.
    expected_visits = 2

    # --- Part 5: Probability of the game ending in finite time for odd n ---
    # For odd n, the random walk of the difference is irreducible on Z_n.
    # With an absorbing state at 0, all other states are transient.
    # The process will be absorbed with probability 1.
    ends_with_prob_one = "yes"
    
    # Print the answers in comma-separated format
    # The problem asks to output the numbers in the equation, so we format the string accordingly.
    
    answer_1 = f"E[X_19] = ({n_19}^2 - 1) / 2 = {int(ex_19)}"
    answer_2 = f"E[X_20] = {str(ex_20)}"
    answer_3 = f"E[X_n] for odd n>1 = {general_formula_odd_n}"
    answer_4 = f"Expected number of times = {expected_visits}"
    answer_5 = f"{ends_with_prob_one}"

    print(f"{answer_1}, {answer_2}, {answer_3}, {answer_4}, {answer_5}")

    # For the final answer format
    final_answer = f"{int(ex_19)}, {str(ex_20)}, {general_formula_odd_n}, {expected_visits}, {ends_with_prob_one}"
    # This part is just for the final answer submission format, not for user output.
    # print(f"<<<{final_answer}>>>")


solve_puzzle()
# The final answer required by the user format is a single comma-separated string.
# E[X_19], E[X_20], E[X_n] for odd n, Expected visits, Finite time?
# 180, inf, (n^2 - 1) / 2, 2, yes