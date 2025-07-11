import math

def solve():
    """
    This function calculates the answers to the user's questions and prints them in a comma-separated format.
    """

    # --- E[X_19] ---
    # For odd n, the formula for the expected time is E[X_n] = (n^2 - 1) / 2.
    n_19 = 19
    e_19 = (n_19**2 - 1) / 2
    # e_19_calc = f"({n_19-1})*({n_19+1})/2 = {int(e_19)}"
    
    # --- E[X_20] ---
    # For even n, the game never ends, so the expected time is infinite.
    e_20 = "inf"
    
    # --- General formula for odd n ---
    # The derived formula is (n-1)(n+1)/2.
    general_formula_e_n = "(n-1)*(n+1)/2"
    
    # --- Expected number of times for odd n > 30 ---
    # This corresponds to the states where the distance between gifts is 11 or n-11.
    # The approximate number of visits is n+1.
    expected_times = "n+1"
    
    # --- Game ends with probability 1 ---
    # For odd n, the absorbing states are reachable from any state, so the game must end.
    ends_prob_one = "yes"
    
    # --- Combine answers ---
    # Give the answers separated by a comma.
    final_answer = f"{int(e_19)},{e_20},{general_formula_e_n},{expected_times},{ends_prob_one}"
    print(final_answer)

solve()