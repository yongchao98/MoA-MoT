def solve_nmr_problem():
    """
    Analyzes the reaction of Pr-DAOTA with sulfuric acid and predicts the
    1H NMR splitting pattern and integration for the most deshielded proton
    in the resulting Compound 1.
    """

    print("--- Analysis of the Reaction and Product ---")
    print("1. The reaction is an electrophilic aromatic sulfonation. Concentrated sulfuric acid adds sulfonic acid (-SO3H) groups to the aromatic rings of the Pr-DAOTA cation.")
    print("2. The introduction of polar -SO3H groups makes the final product, Compound 1, water-soluble as stated in the problem.")
    print("3. The most deshielded proton in Compound 1 is the one on the central ring, located between the two electron-withdrawing nitrogen atoms and within a positively charged system. This environment causes its signal to appear at a very high chemical shift (downfield).")
    print("-" * 40)

    print("--- Determining Integration and Splitting Pattern ---")
    
    # Step 1: Integration
    print("\nIntegration Calculation:")
    integration_value = 1
    print(f"By examining the structure of Compound 1, we can see there is only one proton in this unique chemical environment.")
    print(f"Therefore, the integration for this peak is: {integration_value}H")

    # Step 2: Splitting Pattern
    print("\nSplitting Pattern Calculation (n+1 rule):")
    print("The splitting pattern depends on 'n', the number of protons on adjacent atoms.")
    print("The carbon atom bonded to our proton of interest is connected to two bridgehead carbons in the fused ring system. These bridgehead carbons do not have any protons attached.")
    
    number_of_neighbors = 0
    print(f"Number of neighboring protons (n) = {number_of_neighbors}")
    
    # The 'equation' for the n+1 rule
    splitting_result = number_of_neighbors + 1
    
    pattern = ""
    if splitting_result == 1:
        pattern = "singlet"
    elif splitting_result == 2:
        pattern = "doublet"
    elif splitting_result == 3:
        pattern = "triplet"
    else:
        pattern = "multiplet"
        
    print(f"Using the n+1 rule, the number of peaks in the signal is n + 1 = {number_of_neighbors} + 1 = {splitting_result}.")
    print(f"A signal composed of {splitting_result} peak is called a {pattern}.")
    print("-" * 40)

    print("\n--- Final Answer ---")
    print(f"The highest deshielded proton peak shows a splitting pattern of a {pattern} with an integration of {integration_value}H.")

# Execute the analysis
solve_nmr_problem()