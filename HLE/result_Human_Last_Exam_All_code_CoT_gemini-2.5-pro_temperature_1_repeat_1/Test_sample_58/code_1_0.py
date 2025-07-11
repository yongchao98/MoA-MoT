import sys

def solve_topological_classification():
    """
    Determines the topological invariant group for a defect in a 2D free fermion system.
    """
    # Step 0: Define problem parameters from the user's query
    system_dimension_d = 2
    T_squared = -1
    P_squared = -1
    defect_description = "point defect"
    stated_codimension_D = 1

    print("--- Step 1: Clarifying Problem Parameters ---")
    print(f"The system is a {system_dimension_d}D free fermion model.")
    print(f"The user described the defect as a '{defect_description} (codimension D={stated_codimension_D})'.")
    print("This contains a contradiction: a point defect (0-dimensional) in a 2D system has codimension D = (system dimension) - (defect dimension) = 2 - 0 = 2.")
    print("A defect with codimension D=1 would be a line defect.")
    print("We will assume 'point defect' is the intended feature, as it is more physically specific. Therefore, we will use the corrected codimension D = 2.\n")
    codimension_D = 2

    print("--- Step 2: Identifying the Symmetry Class ---")
    print(f"The system has time-reversal symmetry T with T^2 = {T_squared}.")
    print(f"The system has particle-hole symmetry P with P^2 = {P_squared}.")
    # These symmetries define the Altland-Zirnbauer (AZ) class.
    # T^2=-1, P^2=-1 corresponds to class CII.
    az_class_name = "CII"
    # The real AZ classes are indexed by an integer 's' from 0 to 7. Class CII corresponds to s=5.
    s_index = 5
    print(f"This combination of symmetries corresponds to the AZ class '{az_class_name}', which has the index s = {s_index}.\n")

    print("--- Step 3: Applying the Defect Classification Rule ---")
    print("The classification of a topological defect is given by the bulk classification of a system in a different dimension.")
    print("The rule is: DefectGroup(class s, codimension D) = BulkGroup(class s, dimension d')")
    print("where the effective dimension d' is calculated as: d' = D - 1.")
    
    effective_dimension_d_prime = codimension_D - 1
    
    print(f"For our problem, s = {s_index} and D = {codimension_D}.")
    print(f"The effective dimension is d' = {codimension_D} - 1 = {effective_dimension_d_prime}.\n")

    print("--- Step 4: Finding the Invariant Group from the Periodic Table ---")
    print(f"We need to find the bulk topological classification for class '{az_class_name}' (s={s_index}) in d'={effective_dimension_d_prime} spatial dimension.")
    
    # The periodic table for real topological insulators. Rows are s (0-7), columns are d (0-3).
    # Source: Ryu, Schnyder, Furusaki, Ludwig (2010)
    periodic_table = {
        # s: (d=0, d=1, d=2, d=3)
        0: ('Z', '0', '0', '0'),   # AI
        1: ('Z2', 'Z', '0', '0'),  # BDI
        2: ('Z2', 'Z2', 'Z', '0'), # D
        3: ('0', 'Z2', 'Z2', 'Z'), # DIII
        4: ('Z', '0', 'Z2', 'Z2'), # AII
        5: ('0', 'Z', '0', 'Z2'),  # CII
        6: ('0', '0', 'Z', '0'),   # C
        7: ('0', '0', '0', 'Z')    # CI
    }

    if s_index in periodic_table and 0 <= effective_dimension_d_prime < len(periodic_table[s_index]):
        result_group = periodic_table[s_index][effective_dimension_d_prime]
        print(f"Looking up the entry for s={s_index} and d'={effective_dimension_d_prime} in the periodic table...")
        print(f"The topological invariant is classified by the group: {result_group}")
        # Use a file-like object for the final answer format to avoid printing it with a newline
        sys.stdout.flush()
        answer_io = open(1, 'w', 1) # line-buffered stdout
        answer_io.write(f'<<<{result_group}>>>')
        answer_io.flush()

    else:
        print("Could not determine the classification from the provided table data.")

if __name__ == '__main__':
    solve_topological_classification()