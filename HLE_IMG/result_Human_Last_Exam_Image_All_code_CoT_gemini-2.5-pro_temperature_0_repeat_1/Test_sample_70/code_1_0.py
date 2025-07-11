def describe_reactions():
    """
    Describes the two pericyclic reactions involved in the transformation.
    """
    # Details of the first reaction
    reaction1_electrons = 4
    reaction1_type = "electrocyclic ring opening"
    reaction1_stereochemistry = "conrotatory"

    # Details of the second reaction
    reaction2_electrons = 6
    reaction2_type = "electrocyclic ring closure"
    reaction2_stereochemistry = "disrotatory"

    print("The thermal transformation involves two sequential pericyclic reactions:")
    print(f"1. A {reaction1_electrons}π-electron {reaction1_stereochemistry} {reaction1_type}.")
    print(f"2. A {reaction2_electrons}π-electron {reaction2_stereochemistry} {reaction2_type}.")

# Execute the function to print the answer
describe_reactions()