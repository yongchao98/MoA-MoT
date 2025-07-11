def analyze_proton_nmr():
    """
    Analyzes the NMR properties of the most deshielded proton in Compound 1.
    Compound 1 is formed by the sulfonation of Pr-DAOTA.
    """

    # Step 1 & 2: Identify the proton and its neighbors.
    # The highest deshielded proton is the unique proton on the central, electron-deficient,
    # positively charged aromatic ring.
    # We inspect the structure to find the number of vicinal neighbors (protons on adjacent carbons).
    # The carbons adjacent to the one bearing this proton are part of ring fusions and have no protons.
    num_vicinal_neighbors = 0

    # Step 3: Determine the integration.
    # There is only one such proton in the entire molecule.
    integration_value = 1

    # Step 4: Apply the n+1 rule to determine the splitting pattern.
    multiplicity = num_vicinal_neighbors + 1
    
    # Convert multiplicity number to its common name.
    patterns = {1: "singlet", 2: "doublet", 3: "triplet", 4: "quartet"}
    splitting_pattern = patterns.get(multiplicity, f"{multiplicity}-plet")

    # Print the analysis step-by-step.
    print("Analysis of the highest deshielded proton peak in Compound 1:")
    print(f"1. The proton is on the central cationic ring.")
    print(f"2. Number of vicinal neighbor protons (n) = {num_vicinal_neighbors}")
    print(f"3. Integration = {integration_value}H")
    print(f"4. The splitting pattern is determined by the n+1 rule.")
    print(f"   Splitting multiplicity = n + 1 = {num_vicinal_neighbors} + 1 = {multiplicity}")
    print("\nResult:")
    print(f"The splitting pattern is a {splitting_pattern} and the integration is for {integration_value}H.")

# Execute the analysis function.
analyze_proton_nmr()
