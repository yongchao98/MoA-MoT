def solve_topological_classification():
    """
    Calculates the topological invariant group for a defect in a 2D free fermion system.
    """
    # The periodic table for the 8 real Altland-Zirnbauer classes.
    # The table stores the topological classification (homotopy groups pi_k)
    # for each class 's' in k=0 to k=7 dimensions.
    # pi_k(C_s) gives the classification in k dimensions for class s.
    periodic_table = {
        0: {"name": "AI",   "groups": ["Z", "0", "0", "0", "Z", "Z2", "Z2", "0"]},
        1: {"name": "BDI",  "groups": ["Z2", "Z", "0", "0", "0", "Z", "Z2", "Z2"]},
        2: {"name": "D",    "groups": ["Z2", "Z2", "Z", "0", "0", "0", "Z", "Z2"]},
        3: {"name": "DIII", "groups": ["0", "Z2", "Z2", "Z", "0", "0", "0", "Z"]},
        4: {"name": "AII",  "groups": ["Z", "0", "Z2", "Z2", "Z", "0", "0", "0"]},
        5: {"name": "CII",  "groups": ["0", "Z", "0", "Z2", "Z2", "Z", "0", "0"]},
        6: {"name": "C",    "groups": ["0", "0", "Z", "0", "Z2", "Z2", "Z", "0"]},
        7: {"name": "CI",   "groups": ["0", "0", "0", "Z", "0", "Z2", "Z2", "Z"]},
    }

    # Step 1: Identify the symmetry class from the problem statement.
    # Symmetries T^2=-1 and P^2=-1 uniquely define the class CII.
    # In the 8-fold classification of real classes, CII has index s=5.
    s_index = 5
    class_name = periodic_table[s_index]["name"]
    print("Step 1: Identify Symmetry Class")
    print(f"The system has symmetries T^2 = -1 and P^2 = -1.")
    print(f"This corresponds to the Altland-Zirnbauer class '{class_name}', which has index s = {s_index}.\n")

    # Step 2: Determine the codimension of the defect.
    # The system is 2D, and the defect is a point (0D).
    d_spatial = 2
    d_defect = 0
    codimension = d_spatial - d_defect
    print("Step 2: Determine Defect Codimension")
    print("A point defect is a defect of dimension d_defect = 0.")
    print(f"In a system with spatial dimension d = {d_spatial}, the codimension D is:")
    print(f"D = d - d_defect = {d_spatial} - {d_defect} = {codimension}\n")

    # Step 3: Apply the classification rule for topological defects.
    # The invariant group is given by the homotopy group pi_{D-1}(C_s).
    pi_index = codimension - 1
    print("Step 3: Apply the Classification Rule")
    print("The topological invariant group for a defect is given by the formula: pi_{D-1}(C_s).")
    print(f"Using D = {codimension} and s = {s_index}, we need to find the group pi_{{{codimension}-1}}(C_{s_index}), which is pi_{{{pi_index}}}(C_{s_index}).\n")

    # Step 4: Look up the result in the periodic table.
    result_group = periodic_table[s_index]["groups"][pi_index]
    print("Step 4: Find the Result")
    print(f"Looking up the entry for class {class_name} (row s={s_index}) and homotopy group pi_{pi_index} (column k={pi_index}):")
    print(f"The topological invariant is classified by the group: {result_group}")

solve_topological_classification()
<<<Z>>>