import numpy as np

def solve_utility():
    """
    Solves the game theory problem to find Alice's expected utility.
    """
    # Step 1: Define the payoff matrix.
    # Rows are Alice's actions, Columns are Bob's actions.
    # The values are tuples of (Alice's Payoff, Bob's Payoff).
    # Actions: 0=Rest, 1=Bike, 2=Run
    payoff_matrix = np.array([
        [(0, 0), (2, 0), (4, 0)],   # Alice Rests
        [(0, 2), (-2, -2), (2, 0)],  # Alice Bikes
        [(0, 4), (0, 2), (-3, -3)]   # Alice Runs
    ])

    # Extract Alice's payoffs for her calculations
    alice_payoffs = payoff_matrix[:, :, 0]

    print("This script determines Alice's expected utility based on game theory.")
    print("The core assumptions are that Alice and Bob are superrational and think identically,")
    print("which means they must arrive at the same strategy (a symmetric equilibrium).\n")

    # Step 2: Find and evaluate symmetric pure strategy equilibria.
    # These occur on the diagonal of the matrix.
    diag_payoffs = np.diag(alice_payoffs)
    pure_rest_payoff = diag_payoffs[0]

    print("--- Analyzing Symmetric Pure Strategy Equilibria ---")
    print("Since both players must choose the same action, we look at the diagonal payoffs:")
    print(f"1. If both Rest: Alice's payoff = {diag_payoffs[0]}")
    print(f"2. If both Bike: Alice's payoff = {diag_payoffs[1]}")
    print(f"3. If both Run:  Alice's payoff = {diag_payoffs[2]}\n")

    # The (Rest, Rest) outcome is a Nash Equilibrium, yielding a payoff of 0.
    # The others are not, as players have an incentive to switch.
    # Superrational players can recognize (Rest, Rest) as a stable, symmetric outcome.
    print(f"The only stable pure strategy symmetric equilibrium is for both to Rest.")
    print(f"Utility for Alice in this equilibrium: {pure_rest_payoff}\n")

    # Step 3: Find and evaluate symmetric mixed strategy equilibria.
    # We solve the system of equations for expected utilities.
    # From algebraic analysis, a mixed equilibrium exists between Bike and Run.
    # E[Bike] = p_bike*(-2) + p_run*(2)
    # E[Run]  = p_bike*(0) + p_run*(-3)
    # Setting them equal: -2*p_bike + 2*p_run = -3*p_run => 2*p_bike = 5*p_run
    # With p_bike + p_run = 1, we get p_run = 2/7, p_bike = 5/7.
    p_run = 2/7
    p_bike = 5/7
    
    # Calculate the utility in this mixed equilibrium
    mixed_utility = -3 * p_run

    print("--- Analyzing Symmetric Mixed Strategy Equilibria ---")
    print("A stable mixed strategy is found where Alice and Bob both randomly choose between Bike and Run.")
    print(f"Probabilities: p(Bike) = 5/7, p(Run) = 2/7.\n")
    print("Alice's expected utility in this scenario is calculated from any action in the mix.")
    print("Using the 'Run' action, the final equation is:")
    print(f"Expected Utility = Payoff(Run, Run) * p(Run) + Payoff(Run, Bike) * p(Bike)")
    print(f"Expected Utility = {alice_payoffs[2][2]} * {p_run:.4f} + {alice_payoffs[2][1]} * {p_bike:.4f} = {mixed_utility:.4f}")
    # Or simplified using the equilibrium condition
    #print(f"Expected Utility = Payoff(Run, Run) * p(Run) = -3 * (2/7) = {mixed_utility:.4f}")

    print("\n--- Final Conclusion ---")
    # Step 4: Compare the equilibria.
    # Superrational players choose the equilibrium with the best payoff.
    final_utility = max(pure_rest_payoff, mixed_utility)

    print("Alice and Bob must choose between two possible stable symmetric outcomes:")
    print(f"1. The 'Both Rest' pure strategy, with a utility of {pure_rest_payoff}.")
    print(f"2. The 'Bike/Run Mix' strategy, with an expected utility of {mixed_utility:.4f}.\n")
    print("Being superrational, they will coordinate on the choice with the higher payoff.")
    
    print("The final calculation for Alice's expected utility is:")
    print(f"Alice's Utility = max(Payoff from Pure Equilibrium, Payoff from Mixed Equilibrium)")
    print(f"Alice's Utility = max({pure_rest_payoff}, {mixed_utility:.4f}) = {final_utility}")


solve_utility()