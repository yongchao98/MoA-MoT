def calculate_fst(N, m):
    """
    Calculates the expected Fst at equilibrium in a simple island model.

    Args:
        N (int): The effective population size of the subpopulations.
        m (float): The migration rate (gene flow) between subpopulations.

    Returns:
        float: The calculated Fst value.
    """
    # The term 4Nm represents the number of effective migrants per generation.
    # If gene flow (m) is zero, Fst is 1 (complete differentiation).
    # As gene flow increases, Fst approaches 0 (homogenization).
    if 4 * N * m == -1:
        # Avoid division by zero, though biologically m is >= 0
        return float('inf')
    return 1 / (1 + 4 * N * m)

def main():
    """
    Demonstrates the relationship between gene flow and Fst.
    """
    N = 1000  # Assume an effective population size of 1000

    print("This script demonstrates that Fst (population differentiation) is low when gene flow is high.")
    print(f"Using a constant population size N = {N}.\n")

    # A list of different migration rates (m) to test, from low to high
    migration_rates = [0.0, 0.00001, 0.0001, 0.001, 0.01, 0.1]

    for m in migration_rates:
        # Calculate Nm, the number of effective migrants
        Nm = N * m
        # Calculate Fst
        fst = calculate_fst(N, m)
        
        # Print the equation and result
        print(f"When migration rate m = {m:<8} (gene flow Nm = {Nm:.2f}):")
        # To show the actual numbers in the final equation:
        print(f"  Fst = 1 / (1 + 4 * {N} * {m}) = {fst:.4f}")
        
        # Add an interpretation
        if fst > 0.25:
            print("  Interpretation: Very high differentiation (low gene flow).\n")
        elif fst > 0.15:
            print("  Interpretation: High differentiation.\n")
        elif fst > 0.05:
            print("  Interpretation: Moderate differentiation.\n")
        else:
            print("  Interpretation: Low differentiation (high gene flow).\n")
    
    print("As shown, when gene flow (m) is high, Fst cannot also be high.")

if __name__ == "__main__":
    main()
