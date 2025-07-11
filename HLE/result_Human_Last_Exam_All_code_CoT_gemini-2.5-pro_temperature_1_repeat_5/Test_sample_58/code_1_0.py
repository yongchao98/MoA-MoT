import sys

def solve_topological_classification():
    """
    Calculates the topological invariant group for a specified free fermion system
    based on the tenfold way classification.
    """
    # The tenfold way periodic table is based on K-theory. The classification
    # of defects can be found using this table.

    # 1. Define the Altland-Zirnbauer (AZ) symmetry classes and their integer indices 'q'.
    az_classes = {
        "AI": 0, "BDI": 1, "D": 2, "DIII": 3,
        "AII": 4, "CII": 5, "C": 6, "CI": 7
    }
    # Create a reverse mapping to get class names from their indices.
    az_class_names = {v: k for k, v in az_classes.items()}

    # This is the periodic table for the classification of 0-dimensional systems,
    # which forms the basis for the defect classification. The index of the list
    # corresponds to the class index 'q'. The values are the topological invariant groups.
    # The sequence is Z, Z_2, Z_2, Z_2, 0, 0, Z, 0 for q = 0 to 7.
    pi_0_groups = ["Z", "Z_2", "Z_2", "Z_2", "0", "0", "Z", "0"]

    # --- Problem Parameters ---
    # The problem specifies T^2 = -1 and P^2 = -1. This corresponds to class CII.
    target_class_name = "CII"
    # The system is 2D with a point defect.
    bulk_dimension = 2
    defect_dimension = 0  # A point is 0-dimensional.

    # --- Calculation ---
    print("Step-by-step derivation:")

    # Step 1: Identify the symmetry class 'q'.
    q = az_classes[target_class_name]
    print(f"1. The system has time-reversal and particle-hole symmetries with T^2=-1 and P^2=-1.")
    print(f"   This corresponds to the Altland-Zirnbauer symmetry class '{target_class_name}'.")
    print(f"   The integer index for this class is q = {q}.")
    print("-" * 30)

    # Step 2: Determine the defect codimension 'D'.
    codimension_D = bulk_dimension - defect_dimension
    print(f"2. The system is 2D, and the defect is a point defect (0D).")
    print(f"   The codimension of the defect is D = (bulk dimension) - (defect dimension).")
    print(f"   D = {bulk_dimension} - {defect_dimension} = {codimension_D}.")
    print("   (Note: We interpret 'point defect' as the primary information, meaning D=2, rather than the value D=1 mentioned in the prompt).")
    print("-" * 30)

    # Step 3: Calculate the effective class index 'k'.
    k = (q - codimension_D) % 8
    effective_class_name = az_class_names[k]
    print("3. The topological classification of a codimension-D defect in a system of class 'q'")
    print("   is given by the classification of a 0D system of class 'k', where:")
    print(f"   k = (q - D) mod 8")
    print(f"   k = ({q} - {codimension_D}) mod 8 = {q - codimension_D} mod 8 = {k}")
    print(f"   The class corresponding to the new index k={k} is '{effective_class_name}'.")
    print("-" * 30)

    # Step 4: Look up the final topological group.
    topological_group = pi_0_groups[k]
    print("4. We look up the invariant group for a 0D system of class 'k' in the periodic table.")
    print(f"   The topological invariant group for class {effective_class_name} (k={k}) is: {topological_group}")
    print("=" * 30)

    print(f"Final Answer: The group of the topological invariant is {topological_group}.")

if __name__ == '__main__':
    solve_topological_classification()
