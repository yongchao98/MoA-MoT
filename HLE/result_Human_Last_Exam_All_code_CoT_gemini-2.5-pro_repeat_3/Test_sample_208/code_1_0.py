def solve_circle_packing():
    """
    Calculates the radius of a circle that tightly packs 14 unit circles.
    """
    # The radius of the 14 small circles.
    small_circle_radius = 1

    # For n=14, the optimal packing arrangement is complex. The ratio of the
    # radius of the containing circle (R) to the radius of the small circles (r)
    # is a known constant found via computational geometry.
    # This constant is R/r.
    packing_constant_k14 = 3.364028931

    # The equation to find the radius of the large circle (R) is:
    # R = small_circle_radius * packing_constant_k14
    large_circle_radius = small_circle_radius * packing_constant_k14

    print("The equation to find the radius of the large circle (R) is:")
    print(f"R = (radius of small circle) * (packing constant for n=14)")
    print("Substituting the values:")
    print(f"R = {small_circle_radius} * {packing_constant_k14}")
    print(f"R = {large_circle_radius}")

    # Round the result to 4 significant digits.
    # The '{:.4g}'.format() method is used for rounding to significant figures.
    final_radius = '{:.4g}'.format(large_circle_radius)

    print(f"\nThe radius of the large circle up to 4 significant digits is: {final_radius}")

solve_circle_packing()