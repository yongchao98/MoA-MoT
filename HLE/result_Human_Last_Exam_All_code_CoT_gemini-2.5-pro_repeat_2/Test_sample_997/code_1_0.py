import math

def solve_game():
    """
    This function solves the sequential game for Agent C's optimal strategy.
    """
    # Step 1: Define the game parameters based on the problem description.
    # The game is sequential: C -> B -> A.
    # The winner is the agent with the highest outcome p(w).
    # p(w) is a strictly decreasing function with p(0) = 1.
    # w values must be unique and non-negative.

    # Step 2: Solve using backward induction.

    # Agent A's move (last): A observes p_B and p_C. A wins by choosing p_A > max(p_B, p_C),
    # which is always possible unless max(p_B, p_C) is already 1.

    # Agent B's move (second): B observes p_C. B knows A will win unless B sets an unbeatable
    # maximum. B can do this by choosing p_B = 1 (if p_C < 1).

    # Agent C's move (first): C knows that if C chooses any p_C < 1, B will choose p_B = 1 and win.
    # To win, C must preempt B. C can do this by choosing the highest possible value for p first.
    
    # The highest possible value for p(w) is p(0).
    p_at_zero = 1.0

    # Therefore, Agent C's optimal strategy is to choose w_C = 0.
    # This gives C an unbeatable success probability.
    p_C = p_at_zero
    
    print("Based on backward induction in this sequential game:")
    print("1. Agent A, as the last mover, wins by picking a higher value unless the maximum is already 1.")
    print("2. Agent B, to prevent A from winning, will pick p_B = 1 if possible.")
    print("3. Agent C, as the first mover, preempts B by choosing w_C = 0, which results in the maximum possible success probability.")
    print(f"The maximum success probability is p(0), which is given as {p_at_zero}.")
    print(f"So, Agent C's optimal success probability is p_C = {p_C}")

    # Step 3: Calculate the final result as requested.
    multiplier = 100
    final_value = math.floor(multiplier * p_C)

    print("\nThe problem asks for the floor of 100 * p_C:")
    print(f"floor({multiplier} * {p_C}) = {final_value}")

solve_game()