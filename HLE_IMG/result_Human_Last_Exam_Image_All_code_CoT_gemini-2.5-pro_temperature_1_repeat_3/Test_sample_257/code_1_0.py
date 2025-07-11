def solve_nmr_problem():
    """
    This script analyzes a chemical reaction and predicts the NMR properties of the product.
    
    The reaction is the sulfonation of di-propyl diazaoxatriangulenium (Pr-DAOTA).
    The goal is to find the splitting pattern and integration of the most deshielded proton peak in the product.
    """
    
    # Step 1: Analyze the reaction and product structure
    explanation = [
        "1. The reaction is an electrophilic aromatic substitution (sulfonation) using concentrated sulfuric acid.",
        "2. Symmetrical disulfonation is predicted to occur on the two bottom benzene rings at the positions para to the central oxygen atom. This creates the water-soluble Compound 1.",
        "3. This reaction leaves one proton on each of the two bottom rings. Due to the product's symmetry, these two protons are chemically equivalent.",
        "4. These two protons are the most deshielded in the molecule. This is because they are ortho to three strong electron-withdrawing features: the central cationic carbon, the oxonium oxygen, and the newly introduced sulfonic acid group (-SO3H)."
    ]
    
    print("Chemical Analysis:")
    for step in explanation:
        print(step)
    
    print("\n" + "-"*30 + "\n")
    
    # Step 2: Determine the NMR properties of the most deshielded peak
    integration_value = 2
    splitting_pattern = "singlet"
    
    print("Prediction for the Highest Deshielded Proton Peak in 1H NMR:")
    
    # Step 3: Output the final answer clearly, including the numbers.
    print(f"Integration: The peak corresponds to {integration_value} equivalent protons.")
    print(f"Splitting Pattern: Each proton has no adjacent protons to couple with, resulting in a '{splitting_pattern}'.")
    
    print("\nFinal Answer:")
    print(f"The splitting pattern is a {splitting_pattern} and the integration is {integration_value}H.")

# Execute the function to solve the problem
solve_nmr_problem()
