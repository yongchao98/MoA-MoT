def solve_biology_control_question():
    """
    This script explains the logic for the experimental control in the user's question.
    It breaks down the purpose and timing of the anti-FLAG antibody.
    """
    
    # Experimental parameters
    inhibitor_concentration_mM = 500
    protein_tag = "FLAG"
    
    # Step 1: Define the problem
    print("Problem: We observe that an anti-MUC1-TN antibody shows less binding when 500 mM GalNAc is present.")
    print("Hypothesis: The free GalNAc sugar is competitively inhibiting the antibody.")
    print("-" * 20)

    # Step 2: Identify the potential confounding variable
    print("Potential Issue: A high concentration of a substance like 500 mM GalNAc might be stressful to cells.")
    print("This stress could cause the cells to reduce the amount of MUC1 protein expressed on their surface.")
    print("If this happens, the antibody signal would decrease because its target is gone, not because of competition.")
    print("-" * 20)

    # Step 3: Explain the role of the control
    print(f"The Control: An anti-{protein_tag} antibody.")
    print(f"This antibody binds to the {protein_tag} tag on the MUC1 protein, a part that is separate from the sugar binding site.")
    print("Therefore, its binding directly measures the total amount of MUC1 on the cell surface.")
    print(f"Purpose: To verify that the {inhibitor_concentration_mM} mM GalNAc treatment has not altered the surface expression of MUC1.")
    print("-" * 20)

    # Step 4: Explain the timing
    print("Timing: In an immunoassay, antibodies that bind directly to the target on the cell are called 'primary antibodies'.")
    print("The anti-FLAG antibody is a primary antibody. It must be added during the primary antibody incubation step.")
    print("-" * 20)

    # Step 5: State the final conclusion
    print("Conclusion: The anti-FLAG antibody should be added with the primary antibodies to verify GalNAc has not altered surface expression of MUC1.")
    print("\nBased on this logic, the correct option is C.")

# Execute the reasoning
solve_biology_control_question()
<<<C>>>