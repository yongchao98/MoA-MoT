def solve_nmr_problem():
    """
    This script determines the splitting pattern and integration of the highest
    deshielded proton in Compound 1 based on chemical principles.

    Step 1: Define the properties of the most relevant proton in Compound 1.
    Based on chemical analysis, the reaction is a disulfonation. The most
    deshielded proton is the one on the central position of the top aromatic ring.
    """
    most_deshielded_proton = {
        'description': "Proton on the central position of the top ring",
        'integration': 1,  # There is only one such proton in the molecule.
        'num_neighbors': 2  # It is adjacent to two equivalent protons.
    }

    print("Analysis of the highest deshielded proton peak in Compound 1:")
    print("-" * 55)

    # Step 2: Determine integration.
    integration = most_deshielded_proton['integration']
    print(f"Integration: The peak corresponds to {integration} proton(s), so its integration is {integration}H.")

    # Step 3: Determine the splitting pattern using the n+1 rule.
    n = most_deshielded_proton['num_neighbors']
    multiplicity = n + 1
    
    patterns = {1: 'singlet', 2: 'doublet', 3: 'triplet', 4: 'quartet', 5: 'quintet'}
    splitting_pattern = patterns.get(multiplicity, f'{multiplicity}-plet')

    print("\nSplitting Pattern Calculation (n+1 rule):")
    print(f"The proton has n = {n} equivalent neighboring protons.")
    # The user prompt requested to output each number in the final equation.
    print(f"The multiplicity is calculated as n + 1 = {n} + 1 = {multiplicity}.")
    print(f"A peak with a multiplicity of {multiplicity} is called a '{splitting_pattern}'.")

    print("\n" + "="*20 + " Final Answer " + "="*20)
    print(f"The splitting pattern is a {splitting_pattern}.")
    print(f"The integration is {integration}H.")
    print("="*55)

solve_nmr_problem()