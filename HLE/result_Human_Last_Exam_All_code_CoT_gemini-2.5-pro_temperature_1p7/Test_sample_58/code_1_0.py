import sys

def find_topological_invariant_group():
    """
    Solves for the topological invariant group for a specified fermionic system with a defect.
    """
    
    # Step 1: Define system parameters from the user's problem description.
    system_spatial_dim = 2
    t_squared = -1
    p_squared = -1
    defect_type = "point defect"
    defect_dim = 0  # A point is 0-dimensional

    print("Analysis of a 2D Free Fermion Model with a Point Defect")
    print("="*60)
    print(f"System Properties: {system_spatial_dim}D spatial dimension, free fermion model.")
    print(f"Symmetries: Time-reversal T with T^2 = {t_squared}, Particle-hole P with P^2 = {p_squared}.")
    print(f"Defect: {defect_type} (dimension {defect_dim}).")
    print("-" * 60)

    # Step 2: Identify the Altland-Zirnbauer (AZ) symmetry class.
    # The combination of T and P symmetries with T^2=-1 and P^2=-1 corresponds to class CII.
    az_class = "CII"
    print(f"Step 1: Identify the symmetry class.")
    print(f"The given symmetries place the system in the Altland-Zirnbauer (AZ) class {az_class}.")
    print("\n")

    # Step 3: Calculate the codimension D of the defect.
    # The prompt contains a contradiction: "point defect" implies D=2, but it explicitly states "D=1".
    # We resolve this by using the standard physical definition.
    codimension = system_spatial_dim - defect_dim
    print(f"Step 2: Determine the defect's codimension (D).")
    print(f"A {defect_type} in a {system_spatial_dim}D system has codimension D = (system dimension) - (defect dimension).")
    print(f"Therefore, D = {system_spatial_dim} - {defect_dim} = {codimension}.")
    print("(Note: We are using D=2 based on the 'point defect' description, as the 'D=1' in the prompt is inconsistent.)")
    print("\n")

    # Step 4: Calculate the effective dimension d' for the classification.
    # The rule is that the defect classification in codimension D is given by the bulk
    # classification in dimension d' = D - 1.
    effective_dim = codimension - 1
    print(f"Step 3: Find the effective dimension (d') for classification.")
    print(f"The classification is given by the bulk invariant in dimension d' = D - 1.")
    print(f"Therefore, d' = {codimension} - 1 = {effective_dim}.")
    print("\n")

    # Step 5: Look up the invariant in the periodic table for the identified class and effective dimension.
    # The periodic table for class CII (s=5) has the following invariants for d = 0, 1, 2, ...
    # d':        0    1    2     3     4    5   6   7
    # Class CII: 0,  2Z,   0,   Z2,   Z2,   Z,  0,  0
    periodic_table = {
        "CII": ["0", "2Z", "0", "Z2", "Z2", "Z", "0", "0"]
    }
    
    if az_class in periodic_table:
        # The table is periodic with period 8, so we use modulo operator.
        invariant_group = periodic_table[az_class][effective_dim % 8]
        print(f"Step 4: Look up the topological invariant group.")
        print(f"For class {az_class} in {effective_dim}D, the periodic table gives the group: {invariant_group}.")
    else:
        print(f"Error: Class {az_class} not found in our simplified table.", file=sys.stderr)
        return

    print("="*60)
    print(f"Conclusion: The group of the topological invariant is {invariant_group}.")
    print(f"<<<{invariant_group}>>>")


if __name__ == "__main__":
    find_topological_invariant_group()