def find_topological_invariant_group():
    """
    Determines the topological invariant group for a defect in a 2D free fermion system.
    """
    # Step 1: Define system parameters from the problem description.
    system_dimension = 2
    defect_codimension = 1
    # Symmetries T^2 = -1 and P^2 = -1 define the AZ class CII.
    az_class = "CII"

    print(f"The system is in AZ symmetry class: {az_class} (T^2=-1, P^2=-1).")
    print(f"The spatial dimension of the system is d = {system_dimension}.")
    print(f"The codimension of the defect is D = {defect_codimension}.")
    print("-" * 20)

    # Step 2: Calculate the effective dimension for the defect classification.
    # The classification of a defect with codimension D in d dimensions is
    # given by the classification of a bulk system in d_eff = d - D dimensions.
    effective_dimension = system_dimension - defect_codimension

    print("Calculating the effective dimension for classification:")
    # The final equation as requested.
    print(f"d_eff = d - D = {system_dimension} - {defect_codimension} = {effective_dimension}")
    print("-" * 20)

    # Step 3: Look up the result in the periodic table for class CII.
    # The classification is periodic with a period of 8 in dimension.
    # We store the classification for class CII for dimensions d = 0 to 7.
    # Table for Class CII:
    # d=0 -> 0
    # d=1 -> Z
    # d=2 -> 0
    # d=3 -> Z_2
    # d=4 -> Z_2
    # d=5 -> Z
    # d=6 -> 0
    # d=7 -> 0
    classification_table_CII = {
        0: "0",
        1: "Z",
        2: "0",
        3: "Z_2",
        4: "Z_2",
        5: "Z",
        6: "0",
        7: "0"
    }

    # The lookup dimension is the effective dimension modulo 8.
    lookup_dimension = effective_dimension % 8
    invariant_group = classification_table_CII[lookup_dimension]

    print("Looking up the invariant group in the tenfold way table:")
    print(f"For class {az_class} in effective dimension {effective_dimension}, the topological invariant group is: {invariant_group}")


if __name__ == "__main__":
    find_topological_invariant_group()
