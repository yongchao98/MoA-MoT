def analyze_braveheart_expression():
    """
    Analyzes and explains the expression pattern of the Braveheart lncRNA
    in embryonic stem cells and differentiating heart cells.
    """

    # --- Step 1: Define the entities and their biological roles ---
    cell_state_initial = "Embryonic Stem Cells (ESCs)"
    cell_state_final = "Differentiating Heart Cells"
    gene_role = "A master regulator required to initiate and drive heart development."

    # --- Step 2: Determine expression level in the initial state ---
    # In pluripotent ESCs, a lineage-specific regulator like Braveheart is not yet needed.
    # Therefore, its expression is expected to be very low or off.
    expression_in_escs = "is very low or absent"

    # --- Step 3: Determine the change in expression during differentiation ---
    # To drive the process of becoming a heart cell, the master regulator must be activated.
    # Therefore, its expression must increase.
    expression_change_during_differentiation = "increases"

    # --- Step 4: Print the biological reasoning ---
    print("Step-by-Step Analysis:")
    print(f"1. In the initial state of '{cell_state_initial}', cells are undifferentiated.")
    print(f"   Since Braveheart's role is '{gene_role}', its expression in this initial state {expression_in_escs}.")
    print("-" * 20)
    print(f"2. To transition into '{cell_state_final}', the cardiac development program must be activated.")
    print(f"   This means the expression of the Braveheart regulator {expression_change_during_differentiation} significantly.")
    print("-" * 20)
    print("3. Conclusion: We need an answer that reflects a LOW initial expression and an INCREASE during differentiation.")
    print("\nEvaluating the choices:")
    print("  A. increases in both (Incorrect)")
    print("  B. decreases in both (Incorrect)")
    print("  C. increases in ESCs, decreases in differentiating cells (Incorrect, it's the opposite)")
    print("  D. decreases in embryonic cells and slightly increases in differentiating heart cells")
    print("     - The phrase 'decreases in embryonic cells' is likely a poor way of saying 'is low in ESCs'.")
    print("     - The phrase 'increases in differentiating heart cells' correctly captures the essential change.")
    print("     - This makes it the most plausible choice.")
    print("  E. increases in ESCs, unchanged in differentiating cells (Incorrect)")
    print("  F. does not change (Incorrect)")
    print("-" * 20)
    print("Final Answer based on analysis: D is the best description among the options.")

analyze_braveheart_expression()