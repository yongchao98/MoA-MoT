def solve_sequential_game():
    """
    Solves the sequential game for Agent C's optimal success probability.

    The solution is derived using backward induction, as explained in the steps above.
    The logic demonstrates that due to the sequential nature of the game,
    a pre-emption dynamic exists.
    """

    # Step 1: Analyze Agent A's optimal strategy.
    # Agent A maximizes P(A wins) by choosing the highest available probability p_A.
    # If p=1 is available, A will choose p_A = 1.
    
    # Step 2: Analyze Agent B's optimal strategy.
    # Agent B knows A will choose p_A=1 if possible. To prevent this and win,
    # B must pre-empt A by choosing p_B = 1. This gives B a win probability of 1
    # against A. This holds true as long as C has not chosen p_C = 1.
    
    # Step 3: Analyze Agent C's optimal strategy.
    # Agent C foresees this entire chain of events.
    # If C chooses p_C < 1, B will choose p_B = 1, resulting in P(C wins) = 0.
    # If C chooses p_C = 1, B and A are prevented from choosing 1.
    # This results in P(C wins) = 1.
    # Therefore, C's optimal choice is 1.
    
    p_C = 1.0

    # The problem asks for the value of floor(100 * p_C).
    # The final equation is simply the calculation based on the derived value of p_C.
    result_value = 100 * p_C
    final_answer = int(result_value) # floor is implicit in int() for positive numbers

    print("Agent C's optimal choice for success probability (p_C) is derived through backward induction.")
    print("The logic shows a pre-emptive strategy is optimal for the first mover.")
    print(f"Optimal p_C = {p_C}")
    print("\nFinal calculation:")
    print(f"100 * p_C = 100 * {p_C} = {result_value}")
    print(f"floor(100 * p_C) = floor({result_value}) = {final_answer}")
    
    return final_answer

# Execute the solution
if __name__ == "__main__":
    final_result = solve_sequential_game()
    print("\nFinal Answer:")
    print(f'The value of floor(100 * p_C) is {final_result}.')
    # Final answer format for the platform
    # print(f'<<<{final_result}>>>')