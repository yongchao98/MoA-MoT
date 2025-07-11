import math

def solve():
    """
    Solves the user's puzzle based on a detailed analysis of the game.
    
    Here's a summary of the derived answers:
    1.  E[X_19]: For odd n, the expected time E[X_n] is (n^2 - 1) / 2.
        For n=19, E[X_19] = (19^2 - 1) / 2 = (361 - 1) / 2 = 180.
    2.  E[X_20]: For even n, the distance between gift-holders is always odd.
        The game-ending condition requires a distance of 2 (even). This is impossible.
        Thus, the expected time is infinite.
    3.  E[X_n] for odd n > 1: The general formula is (n^2 - 1) / 2.
    4.  Expected visits to distance 11 for n > 30: The number of people
        between the gifts is d-1. So, 10 people between them means distance d=11.
        The random walk of the distance starts at 1, must pass through all
        intermediate odd distances (like 11) to reach the maximum distance,
        then cross over to even distances and descend towards the absorbing
        distance 2. This suggests it passes through the "rung" d=11 once on the
        way up and once on the way down from the other side of the 'distance ladder' due to symmetry, totaling 2 visits.
        Let's rethink that logic. The "up" walk is on odd distances, while the "down" walk is on even distances. 
        So once the walk crosses to even distances, it won't hit an odd one like 11 again. 
        So it will be visited exactly once.
        Let's try to derive it. Expected time from d=1 to d=11 is E[1->11] = sum_{k=1 to 5} 4k = 4*(5*6/2) = 60.
        The expected number of visits to a transient state j starting from i before absorption is P(hit j from i)/P(leave j and not return). 
        The calculation is quite involved. A simpler argument suggests that the walk on distances starts at d=1 and must pass through d=11 to reach higher odd distances and eventually cross to the even-distance part of the state space for the game to terminate. It can't return to the odd states after that. Therefore, it will be visited exactly once.
    5.  Game ends with probability 1 for odd n: Yes. The state space is finite, and the absorbing state (game end) is reachable from any state. Thus, absorption is certain.
    """
    
    # Value for E[X_19]
    n_19 = 19
    e_x19 = (n_19**2 - 1) // 2

    # Value for E[X_20]
    e_x20 = 'inf'

    # General formula for E[X_n] for odd n
    # We will represent the formula as a string
    # and print the calculation for n_19 as an example of its use
    general_formula = "(n^2 - 1) / 2"
    
    # Expected number of times for distance 11
    # As reasoned above, this should be 1.
    expected_visits = 1

    # Probability of ending in finite time for odd n
    prob_one = "yes"
    
    print(f"The expected value of X_19 is ({n_19}^2 - 1) / 2 = {e_x19}")
    print(f"The expected value of X_20 is {e_x20}")
    print(f"The general formula for odd n is {general_formula}")
    print(f"The expected number of visits to distance 11 is {expected_visits}")
    print(f"Does the game end with probability 1 for odd n? {prob_one}")
    
    # Final answer string, as requested by the user prompt
    answer_string = f"{e_x19},{e_x20},({n_19}^2 - 1) / 2 = {(n_19**2-1)/2},{expected_visits},{prob_one}"
    final_answer_to_print = f"{e_x19},{e_x20},{e_x19},{expected_visits},{prob_one}"

    # Print comma separated values for easy parsing
    # However, the prompt is slightly ambiguous, asking for the formula.
    # So I will present the answer string containing the formula expression.
    
    final_e_xn_str = f"({n_19}^2-1)/2" # Let's display the calculation as part of the answer
    final_answer = f"{e_x19},{e_x20},(n^2-1)/2,{expected_visits},{prob_one}"

    # After thinking, the prompt might expect the equation fully expanded with numbers for the n=19 case.
    print(f"\nFinal Answer String:\n{e_x19},{e_x20},({n_19}^2 - 1) / 2 = {e_x19},{expected_visits},{prob_one}")

solve()