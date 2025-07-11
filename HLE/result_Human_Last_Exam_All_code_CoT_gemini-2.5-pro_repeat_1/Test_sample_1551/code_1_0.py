import math

def solve():
    """
    Solves the user's puzzle based on the derived formulas.
    """

    # --- Question 1: E[X_19] ---
    # For odd n, E[X_n] = (n^2 - 1) / 2
    n_19 = 19
    ex_19 = (n_19**2 - 1) / 2
    
    # --- Question 2: E[X_20] ---
    # For even n, the game never ends.
    n_20 = 20
    ex_20 = math.inf

    # --- Question 3: E[X_n] for odd n ---
    # The formula is (n^2 - 1) / 2. We'll represent this as a string.
    ex_n_odd_formula = "(n^2 - 1) / 2"

    # --- Question 4: Expected visits for n > 30 ---
    # The distance between gifts is 10, so the number of edges is 11.
    # This corresponds to D_t = 11 or D_t = n - 11.
    # The expected number of visits is 2 * (n - 2*k + 1) where k=11.
    # Let's represent this as a function of n.
    k = 11
    expected_visits_formula = f"2 * (n - 2*{k} + 1) = 2 * (n - 21)"
    # The final equation is 2n - 42. So we output 2, n, -42
    
    # --- Question 5: Does the game end with probability 1 for odd n? ---
    # Yes, because the difference walk D_t is irreducible on a finite state space.
    ends_prob_1 = "Yes"

    # --- Print the results ---
    # For question 3, we show the formula.
    # For question 4, we show the final equation parts.
    print(f"E[X_19]: {int(ex_19)}")
    print(f"E[X_20]: {ex_20}")
    print(f"E[X_n] for odd n: {ex_n_odd_formula}")
    print(f"Expected number of times distance is 10 for odd n > 30: The final equation is 2*n - 42. The numbers are 2, 1, -42.")
    print(f"Game ends with probability 1 for odd n: {ends_prob_1}")

    # The final answer in the requested format
    # The formula for Q4 is 2n - 42. The numbers are 2, -42. But the question asks for 10 and n-12...
    # Let me re-read Q4. "between the two gifts there are exactly 10 and n-12 friends"
    # This means the number of edges are 11 and n-11. This is just one condition, described in two ways.
    # The problem asks to output each number in the final equation. 2*n-42
    final_answer = f"<<<{int(ex_19)}, {ex_20}, (n^2-1)/2, 2(n-21), {ends_prob_1}>>>"
    # I should just output the values separated by comma.
    # For Q4, I'll output the formula string as requested by the user format.
    final_answer_string = f"{int(ex_19)},{ex_20},(n^2-1)/2,2(n-21),{ends_prob_1}"
    print("\nFinal Answer:")
    print(final_answer_string)


solve()