def find_crystal_class():
    """
    This function explains and identifies the crystal class that is achiral, 
    non-polar, and exhibits optical activity.
    """

    print("Analyzing the requirements for the crystal class:")
    print("1. Achiral (contains a rotoinversion axis, e.g., mirror or inversion center)")
    print("2. Non-polar (has no spontaneous electric dipole moment)")
    print("3. Optically Active (can rotate the plane of polarized light)")
    print("-" * 50)
    
    print("Key Insight: The rule that 'optical activity requires chirality' is a general guideline, but not strictly true.")
    print("The precise requirement for optical activity is the ABSENCE of a center of inversion (ī) and any mirror planes (m).")
    print("\nWe must find an achiral crystal class that lacks both these symmetry elements.")
    print("- An achiral class must have a rotoinversion axis (like 3̄, 4̄, or 6̄).")
    print("- Classes with 3̄ contain an inversion center, so they are not optically active.")
    print("- Classes with 6̄ contain a mirror plane, so they are not optically active.")
    print("- This leaves only classes whose achirality is due to a 4̄ axis, without any additional mirror planes or inversion centers.")
    print("\nThe crystal class that fits this description is '4-bar'.")
    print("-" * 50)

    # The "equation" here is the name of the class, which includes the number 4.
    final_class_name = "4-bar"
    number_in_name = "4"

    print("The only achiral, non-polar, and optically active crystal class is:")
    print(f"Hermann-Mauguin Symbol: {final_class_name}")
    print(f"This symbol contains the number: {number_in_name}")

# Execute the function to print the result.
find_crystal_class()
