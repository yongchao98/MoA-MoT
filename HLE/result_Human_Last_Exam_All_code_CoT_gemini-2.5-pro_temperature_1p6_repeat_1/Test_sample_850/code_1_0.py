def find_crystal_classes():
    """
    Analyzes the symmetry requirements for optical activity and achirality
    to answer the user's question.
    """
    
    # Condition 1: Optical Activity
    # A crystal must be chiral to be optically active.
    # Chirality means the crystal's point group LACKS a center of inversion (i)
    # and mirror planes (m).
    condition_for_optical_activity = "must be chiral (lacks i and m)"

    # Condition 2: Achirality
    # An achiral crystal is superimposable on its mirror image.
    # This means its point group MUST HAVE a center of inversion (i) or a
    # mirror plane (m).
    condition_for_achirality = "must be achiral (has i or m)"
    
    print("Analyzing the requirements for crystal classes...")
    print("-" * 50)
    print(f"Requirement for Optical Activity: The crystal class {condition_for_optical_activity}.")
    print(f"Requirement for being Achiral: The crystal class {condition_for_achirality}.")
    print("-" * 50)
    print("\nConclusion:")
    print("The conditions for being 'optically active' and 'achiral' are mutually exclusive.")
    print("A crystal class cannot simultaneously lack and possess an inversion center or mirror plane.")
    print("\nTherefore, there are no crystal classes that are both achiral and optically active.")
    
# Execute the function to print the explanation.
find_crystal_classes()
