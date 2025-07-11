def solve_game():
    """
    This function solves the sequential game theory problem based on backward induction.

    The game is defined as follows:
    1. Three agents (C, B, A) choose a non-negative, unique real number `w` sequentially.
    2. The agent with the smallest `w` wins. This is because their `p(w)` is the highest,
       and `p(w)` is a strictly decreasing function.
    3. The objective for each agent is to maximize their probability of winning.

    The solution is found using backward induction:
    - Agent A's analysis: A moves last. A sees w_B and w_C. If min(w_B, w_C) > 0,
      A can choose w_A < min(w_B, w_C) and win. If min(w_B, w_C) = 0, A cannot win.
    - Agent B's analysis: B sees w_C. B knows A's strategy. If C chooses w_C > 0, B can
      choose w_B = 0. This makes it impossible for A to win, and since w_B < w_C, B wins.
      If C chooses w_C = 0, B cannot win.
    - Agent C's analysis: C moves first. C knows B's and A's strategies. If C chooses
      w_C > 0, B will choose w_B = 0 and win. C will lose. If C chooses w_C = 0,
      B cannot choose 0, so B must choose w_B > 0. A also cannot win. C will win.

    Therefore, C's optimal strategy is to choose w_C = 0.
    """

    # According to the game analysis, C's optimal strategy is to choose w_C = 0.
    w_C = 0

    # The problem defines p(w) such that p(0) = 1.
    # So, C's chosen success probability is p_C.
    p_C = 1

    # The problem asks for the floor of 100 * p_C.
    result = int(100 * p_C)

    print(f"Step 1: Agent C, moving first, analyzes the game.")
    print(f"Step 2: C realizes that if they choose any w_C > 0, Agent B can choose w_B = 0 and win.")
    print(f"Step 3: To prevent this, C's optimal strategy is to preemptively choose the unbeatable value, w_C = {w_C}.")
    print(f"Step 4: This choice corresponds to a success probability p_C = p({w_C}) = {p_C}.")
    print(f"Step 5: The final calculation is floor(100 * p_C).")
    print(f"floor(100 * {p_C}) = {result}")

solve_game()