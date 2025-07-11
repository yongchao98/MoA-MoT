def solve_poker_problem():
    """
    This function formats and prints the solution to the poker problem based on the analysis.
    The logic is as follows:
    - With QQ (Case B), we are beaten and any bet will be called. We must check to minimize loss.
    - With AA (Case A), we have the better hand. However, the villain knows we would only check with QQ.
      Therefore, any bet we make signals that we have AA. A perfect villain will fold to any bet.
      Since both checking and betting win the same amount (the pot), we are indifferent to our action.
      We can choose any strategy. We will provide a mixed strategy as an example.
    - Villain's Response (Case C): Since any bet comes from AA, the villain (with KK) will always fold.
    """
    
    # Define the strategies for each case. All percentages and sizings are even numbers.
    # The mixed strategy for AA is sorted by the amount of money put in ($0 for CHECK, $10 for BET).
    answer_a = "CHECK 50%, BET $10 50%"
    answer_b = "CHECK 100%"
    answer_c = "FOLD 100%"

    # Format the final output string.
    final_output = f"A) {answer_a}\nB) {answer_b}\nC) {answer_c}"
    
    print(final_output)

solve_poker_problem()