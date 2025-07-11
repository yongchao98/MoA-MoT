def find_topological_invariant_group():
    """
    Calculates the topological invariant group for a 2D free fermion system
    with a point defect and specific symmetries.
    """
    # Step 1: Define system parameters from the problem description.
    # We are considering a point defect (0-dimensional) in a 2D system.
    system_dimension = 2
    defect_dimension = 0  # For a point defect
    
    print(f"--- Analyzing a point defect in a {system_dimension}D system ---")
    print("\nStep 1: Determine the codimension of the defect.")
    print(f"The dimension of the bulk system is d_bulk = {system_dimension}.")
    print(f"The dimension of a point defect is d_defect = {defect_dimension}.")
    
    codimension = system_dimension - defect_dimension
    
    print(f"The codimension 'D' is calculated as: D = d_bulk - d_defect = {system_dimension} - {defect_dimension} = {codimension}.")
    print("Note: We are using D=2 based on the physical description 'point defect in 2D', as this is standard.")
    
    # Step 2: Determine the Altland-Zirnbauer (AZ) symmetry class.
    T_squared = -1
    P_squared = -1
    
    print("\nStep 2: Determine the symmetry class.")
    print(f"The system has time-reversal symmetry with T^2 = {T_squared}.")
    print(f"The system has particle-hole symmetry with P^2 = {P_squared}.")

    # For fermions, chiral symmetry S=TP squares to S^2 = -T^2 * P^2
    S_squared = -T_squared * P_squared
    print(f"The corresponding chiral symmetry squares to S^2 = -(T^2) * (P^2) = -({T_squared}) * ({P_squared}) = {S_squared}.")
    
    # Determine the AZ class from the given symmetries.
    # The combination T^2=-1, P^2=-1 corresponds to class CII.
    az_class = "CII"
    print(f"This set of symmetries belongs to the Altland-Zirnbauer class '{az_class}'.")

    # Step 3: Find the classification using the periodic table.
    # The classification of a defect of codimension D is the same as the
    # classification of a bulk system of dimension d_eff = D.
    effective_dimension = codimension
    
    print(f"\nStep 3: Find the topological invariant group.")
    print(f"The classification for a defect of codimension D={codimension} is equivalent to the classification for a bulk system of dimension d_eff = {effective_dimension}.")
    
    # Periodic table for topological insulators (real classes).
    # We only need the entry for class CII at d=2.
    # Source: C.-K. Chiu et al., Rev. Mod. Phys. 88, 035005 (2016).
    periodic_table = {
        "CII": {
            0: "0",      # d=0
            1: "Z",      # d=1
            2: "Z_2",    # d=2
            3: "Z_2"     # d=3
        }
    }
    
    # Look up the result in the table.
    if az_class in periodic_table and effective_dimension in periodic_table[az_class]:
        result_group = periodic_table[az_class][effective_dimension]
    else:
        result_group = "Unknown"

    print(f"\n--- Final Result ---")
    print(f"According to the periodic table, the classification for class '{az_class}' in dimension d={effective_dimension} is the group {result_group}.")
    print(f"Therefore, the group of its topological invariant is {result_group}.")

if __name__ == '__main__':
    find_topological_invariant_group()