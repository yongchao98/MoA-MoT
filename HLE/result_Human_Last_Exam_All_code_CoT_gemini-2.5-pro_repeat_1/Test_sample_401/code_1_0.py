def solve_composant_problem():
    """
    Solves the topological problem regarding the minimum number of composants
    in an indecomposable continuum.

    This script explains the reasoning and prints the final answer.
    """

    # --- Step 1: Define the terms ---
    print("Problem: What is the smallest number of composants an indecomposable continuum can have?")
    print("-" * 80)
    print("1. An indecomposable continuum is a compact, connected space that cannot be split into two smaller, proper subcontinua.")
    print("2. A composant is a specific type of dense sub-set within an indecomposable continuum. The collection of all composants partitions the space.")
    print("-" * 80)

    # --- Step 2: Establish a lower bound ---
    print("Reasoning for the lower bound:")
    print("Can the number be 1?")
    print("No. By definition, a composant of a space X is a union of *proper* subcontinua of X.")
    print("A proper subcontinuum is strictly smaller than X.")
    print("If there were only one composant, it would have to be the entire space X.")
    print("This leads to a contradiction: X would have to be a proper subset of itself.")
    print("Therefore, the number of composants must be at least 2.")
    print("-" * 80)

    # --- Step 3: State the known minimum ---
    print("Finding the achievable minimum:")
    print("For standard 'metric' continua, the number of composants is uncountably infinite.")
    print("However, the problem allows for 'non-metric' continua.")
    print("It is a famous result in topology (by M.E. Rudin, assuming the Continuum Hypothesis) that an indecomposable continuum with exactly 2 composants can be constructed.")
    print("-" * 80)

    # --- Step 4: Final Answer ---
    smallest_number = 2
    print("Conclusion: The smallest possible number of composants is 2.")
    print(f"Final Answer Equation: Smallest number of composants = {smallest_number}")


if __name__ == "__main__":
    solve_composant_problem()
