import sys

def solve_fermion_classification():
    """
    Solves for the topological invariant group of a 2D free fermion system
    with specified symmetries and a defect.
    """
    # Step 1: Define the system parameters from the problem description.
    spatial_dimension = 2
    T_squared = -1
    P_squared = -1
    codimension = 1

    # Step 2: Determine the Altland-Zirnbauer (AZ) symmetry class.
    # The combination T^2=-1 and P^2=-1 corresponds to class CII.
    # We can use a dictionary to map symmetry properties to class names.
    az_class_map = {
        (-1, -1): "CII"
        # Other classes could be added here for a more general tool.
    }
    az_class = az_class_map.get((T_squared, P_squared), "Unknown")

    if az_class == "Unknown":
        print("Could not determine the AZ class for the given symmetries.", file=sys.stderr)
        return

    # Step 3: Use the periodic table for topological classification.
    # The classification of a defect of codimension D is the same as the
    # classification of a bulk system of dimension d=D in the same class.
    # We store the table as a dictionary where keys are AZ classes and
    # values are lists of invariants for d=0 to d=7.
    # Z = integers, Z2 = integers modulo 2, 0 = trivial.
    periodic_table = {
        "CII": ["0", "Z", "0", "Z2", "Z2", "Z", "0", "0"]
    }

    # The effective dimension for the defect classification is its codimension.
    effective_dimension = codimension

    # Look up the invariant group in the table.
    # The list is 0-indexed, so the index matches the dimension.
    if az_class in periodic_table and 0 <= effective_dimension < len(periodic_table[az_class]):
        invariant_group = periodic_table[az_class][effective_dimension]
    else:
        print(f"Classification for class {az_class} at d={effective_dimension} is not available.", file=sys.stderr)
        return

    # Step 4: Print the detailed analysis and the final result.
    print("Analysis of the Fermion System:")
    print("=" * 35)
    print(f"1. System Properties:")
    print(f"   - Bulk Spatial Dimension (d): {spatial_dimension}")
    print(f"   - Time-Reversal Symmetry (T^2): {T_squared}")
    print(f"   - Particle-Hole Symmetry (P^2): {P_squared}")
    print(f"   - Defect Codimension (D): {codimension}")
    print("-" * 35)
    print(f"2. Symmetry Classification:")
    print(f"   The combination of T^2 = {T_squared} and P^2 = {P_squared} places the system in AZ class '{az_class}'.")
    print("-" * 35)
    print(f"3. Defect Classification Principle:")
    print(f"   The topological invariant for a defect of codimension D = {codimension} is given by the")
    print(f"   classification of a bulk system with dimension d_eff = D = {effective_dimension} in the same class '{az_class}'.")
    print("-" * 35)
    print(f"4. Final Result from Periodic Table:")
    print(f"   Looking up the entry for class '{az_class}' and dimension d_eff = {effective_dimension}:")
    print(f"   The group of the topological invariant is: {invariant_group}")
    print("=" * 35)

    # Final answer in the required format
    print(f"\n<<<{invariant_group}>>>")

# Execute the function
solve_fermion_classification()