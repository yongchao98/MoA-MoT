def explain_fuzzy_dimensions():
    """
    This function explains and demonstrates the dimensional structure
    of Type-1, Type-2, and Type-3 fuzzy membership functions.
    """
    x_input = 7.0

    print("--- Illustrating the Dimensional Structure of Fuzzy Membership Functions ---")
    print(f"Let's consider a sample input value x = {x_input}\n")

    # --- Type-1 Fuzzy Set ---
    # Membership is a single number (a point).
    # This is a 2D relationship: (input) -> (membership_grade)
    t1_membership = 0.75
    print("1. Type-1 Membership Function (2-Dimensional Structure)")
    print(f"   - For input x = {x_input}, the membership is a single crisp value: {t1_membership}")
    print("   - This defines a point in a 2D space (Input Axis, Membership Axis).\n")


    # --- Type-2 Fuzzy Set ---
    # Membership is a 1D interval, a slice of the 2D "Footprint of Uncertainty" (FOU).
    # The overall model requires 3 dimensions:
    # (1) input x, (2) primary membership u, (3) secondary membership grade.
    t2_membership_footprint_slice = (0.6, 0.9)
    print("2. Type-2 Membership Function (3-Dimensional Structure)")
    print(f"   - For input x = {x_input}, the membership is an interval (e.g., from {t2_membership_footprint_slice[0]} to {t2_membership_footprint_slice[1]}).")
    print("   - This interval represents the uncertainty and is a 1D slice of a 2D 'Footprint of Uncertainty' (FOU).")
    print("   - The entire model is 3D, describing a surface of uncertainty over the input domain.\n")


    # --- Type-3 Fuzzy Set ---
    # Membership is a Type-2 Set, which can be visualized as a 2D area for a given x.
    # This 2D area is a slice of a 3D "Volume of Uncertainty" (VOU).
    # The overall model requires 4 dimensions.
    # The key change from Type-2 is that the uncertainty model itself becomes 3D.
    t3_membership_slice_example = "a Type-2 Fuzzy Set (a 2D footprint)"
    print("3. Type-3 Membership Function (4-Dimensional Structure)")
    print(f"   - For input x = {x_input}, the membership is a complete {t3_membership_slice_example}.")
    print("   - This creates a 'Volume of Uncertainty' (VOU), which is a 3D entity.")
    print("   - The fundamental difference from Type-2 is this move from a 2D FOU to a 3D VOU.")
    print("   - Therefore, a three-dimensional uncertainty model is added to the structure.\n")


# Run the explanation
explain_fuzzy_dimensions()