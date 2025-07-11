import math

def solve_simulation_cost():
    """
    Calculates the minimal resource cost to simulate the non-local
    correlations of a singlet state using a PR-Box.

    The method is to find the minimum fraction of a PR-Box, when mixed with a
    local classical strategy, required to reproduce the maximal CHSH
    violation achieved by a quantum singlet state.
    """

    # 1. Define the maximal CHSH inequality values for the different models.
    # The CHSH value is a measure of non-locality in the correlations.
    S_lhv = 2.0  # Maximal value for any Local Hidden Variable (LHV) model.
    S_quantum = 2 * math.sqrt(2)  # Maximal value for a quantum singlet state (Tsirelson's bound).
    S_pr_box = 4.0  # Value for a Popescu-Rohrlich (PR) Box, the maximally non-local resource.

    # 2. Set up the equation for the mixture.
    # We simulate the quantum result (S_quantum) with a convex combination of
    # a local model (S_lhv) and a PR-box (S_pr_box). 'p' is the fraction
    # of the PR-box resource required.
    # Equation: S_quantum = (1 - p) * S_lhv + p * S_pr_box

    # 3. Solve for 'p'.
    # p * S_pr_box - p * S_lhv = S_quantum - S_lhv
    # p * (S_pr_box - S_lhv) = S_quantum - S_lhv
    # p = (S_quantum - S_lhv) / (S_pr_box - S_lhv)
    p = (S_quantum - S_lhv) / (S_pr_box - S_lhv)

    # 4. Print the result and the equation.
    # The output explains the context and shows the final equation with all numbers,
    # as requested.
    print("The minimal resource cost can be quantified by the fraction of a PR-Box ('p')")
    print("needed to simulate the non-locality of a quantum singlet state.")
    print("This is calculated by mixing a local model with a PR-Box to match the")
    print("quantum CHSH value.")
    print("\nThe equation for the required PR-Box fraction 'p' is:")
    print(f"p = (S_quantum - S_lhv) / (S_pr_box - S_lhv)")
    print("\nPlugging in the values:")
    print(f"p = ({S_quantum:.3f} - {S_lhv:.3f}) / ({S_pr_box:.3f} - {S_lhv:.3f})")
    print(f"p = {S_quantum - S_lhv:.3f} / {S_pr_box - S_lhv:.3f}")
    print(f"p = {p:.3f}")
    print("\nThis means that a singlet state contains about 41.4% of the non-locality of a single PR-Box.")
    print("\nIt's also known that 1 bit of communication is sufficient to perfectly simulate all correlations")
    print("of the singlet state. The value 'p' above can be interpreted as the average communication")
    print("cost (in bits) of a mixed strategy that simulates the CHSH violation, making it the minimal resource cost under this interpretation.")

if __name__ == '__main__':
    solve_simulation_cost()
    # The final numerical answer is p = sqrt(2) - 1
    final_answer = math.sqrt(2) - 1
    # print(f"\n<<<{final_answer:.3f}>>>") # For extraction