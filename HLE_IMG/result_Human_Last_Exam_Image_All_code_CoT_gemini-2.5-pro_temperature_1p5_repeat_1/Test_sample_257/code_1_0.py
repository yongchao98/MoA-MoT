def solve_chemistry_nmr_problem():
    """
    Analyzes the sulfonation of Pr-DAOTA and determines the 1H NMR
    characteristics of the most deshielded proton in the product.
    """

    print("Step-by-step analysis of the NMR properties of Compound 1:")
    print("-" * 60)

    # Step 1: Identify the most deshielded proton in the product (Compound 1).
    print("1. Identifying the Target Proton:")
    print("   The reaction is an electrophilic sulfonation on the most activated aromatic ring,")
    print("   which is the one with the two N-propyl groups. Symmetrical disulfonation occurs.")
    print("   The most deshielded proton is the single proton remaining on this sulfonated ring,")
    print("   due to its proximity to two N atoms and two electron-withdrawing -SO3H groups in a cationic system.")
    print("")

    # Step 2: Determine the integration of the peak for this proton.
    print("2. Determining the Integration:")
    integration = 1
    print(f"   The product molecule is symmetrical, and there is only one proton in this specific")
    print(f"   chemical environment. Therefore, the integration value is {integration}H.")
    print("")

    # Step 3: Determine the splitting pattern for this proton using the n+1 rule.
    print("3. Determining the Splitting Pattern:")
    print("   The splitting pattern is determined by the 'n+1 rule', where 'n' is the number")
    print("   of protons on adjacent atoms.")
    
    n = 0
    print(f"   For the proton in question, the adjacent carbon atoms are part of ring fusions and have")
    print(f"   no attached protons. Thus, n = {n}.")

    # Apply the n+1 rule calculation
    num_peaks = n + 1
    pattern = "singlet"
    print(f"\n   The number of peaks in the signal is calculated by the equation: n + 1")
    print(f"   Substituting the value of n: {n} + 1 = {num_peaks}")
    print(f"   A signal with {num_peaks} peak is called a {pattern}.")
    print("-" * 60)

    # Final Conclusion
    print("\nFinal Answer:")
    print(f"The splitting pattern of the highest deshielded proton peak is a {pattern},")
    print(f"and its integration is {integration}H.")

# Run the analysis function
solve_chemistry_nmr_problem()