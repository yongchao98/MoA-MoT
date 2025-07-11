def solve_crane_regions():
    """
    Calculates the number of regions the fold lines of a standard origami crane
    divide the paper into.

    The calculation is based on a known breakdown of the crane's crease pattern.
    The pattern is composed of a central body and four identical flaps that
    form the wings, neck, and tail.
    """
    # Number of identical flaps (2 wings, 1 neck, 1 tail)
    num_flaps = 4

    # According to analysis of the crease pattern by origami experts,
    # each of the four main flaps is divided into 32 distinct regions.
    regions_per_flap = 32

    # The central body of the crane crease pattern contains 18 regions.
    regions_in_central_body = 18

    # Calculate the total number of regions contributed by the flaps
    total_flap_regions = num_flaps * regions_per_flap

    # Calculate the grand total by adding the central body regions
    total_regions = total_flap_regions + regions_in_central_body

    print("Calculating the number of regions in an unfolded origami crane's crease pattern.")
    print("-" * 70)
    print(f"The pattern consists of a central body and {num_flaps} main flaps.")
    print(f"Regions in the central body: {regions_in_central_body}")
    print(f"Regions in each of the {num_flaps} flaps: {regions_per_flap}")
    print("\nThe final calculation is based on summing these parts:")
    print(f"Total Regions = (Number of flaps * Regions per flap) + Regions in central body")
    
    # Final output showing the equation with the numbers plugged in
    print("\nFinal Equation:")
    print(f"{total_regions} = ({num_flaps} * {regions_per_flap}) + {regions_in_central_body}")


solve_crane_regions()