def calculate_synthesis_yield():
    """
    Calculates and compares the theoretical yields for direct vs. fragmented
    peptide synthesis to illustrate the benefit of Native Chemical Ligation (NCL).
    """
    # --- Parameters ---
    # A very optimistic coupling efficiency for each amino acid addition.
    coupling_efficiency = 0.99
    # Length of the full peptide.
    full_length = 100
    # Length of a representative fragment for NCL.
    fragment_length = 50

    # --- Calculations ---
    # The number of coupling reactions is the length - 1.
    couplings_full = full_length - 1
    couplings_fragment = fragment_length - 1

    # Calculate the overall theoretical synthesis yield.
    yield_full = coupling_efficiency ** couplings_full
    yield_fragment = coupling_efficiency ** couplings_fragment

    # --- Output ---
    print("Comparing Theoretical Yields of Peptide Synthesis Strategies:\n")
    print(f"Assuming a per-step coupling efficiency of {coupling_efficiency * 100}%\n")

    print("1. Direct Synthesis (for a 100-amino-acid peptide):")
    # The prompt requires showing the numbers in the equation.
    print(f"   Equation: {coupling_efficiency}^{couplings_full}")
    print(f"   Theoretical Yield: {yield_full:.4f} (or {yield_full:.2%})")
    print("-" * 30)

    print("2. Fragment Synthesis (for a 50-amino-acid fragment used in NCL):")
    print(f"   Equation: {coupling_efficiency}^{couplings_fragment}")
    print(f"   Theoretical Yield: {yield_fragment:.4f} (or {yield_fragment:.2%})")
    print("-" * 30)

    print("\nConclusion:")
    print(f"The synthesis yield for a 50aa fragment ({yield_fragment:.2%}) is significantly higher")
    print(f"than for a 100aa peptide ({yield_full:.2%}). This makes synthesizing two smaller fragments")
    print("and joining them via NCL a much more practical and higher-yielding approach.")

calculate_synthesis_yield()