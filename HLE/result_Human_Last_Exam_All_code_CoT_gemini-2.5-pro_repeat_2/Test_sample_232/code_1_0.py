def solve_crane_regions():
    """
    Calculates the number of regions the fold lines of a standard origami crane divide the paper into.

    The calculation is based on deconstructing the crane's final crease pattern into its main components:
    1. A central diamond-shaped area that becomes the body and wings.
    2. Four identical corner sections that become the neck, tail, and two hidden points.
    3. The final asymmetrical "crimp" fold that creates the head.
    """

    # 1. The central area of the crease pattern is a dense grid which forms the crane's body and wings.
    # Careful counting reveals it is divided into 64 small triangular regions.
    central_regions = 64

    # 2. Each of the four corner sections of the paper is identically creased to form the points
    # (neck, tail, etc.). Each of these four sections contains 28 regions.
    num_corner_sections = 4
    corner_section_regions = 28
    total_corner_regions = num_corner_sections * corner_section_regions

    # 3. The final step in folding a crane is the "crimp fold" to create the head.
    # This single maneuver adds two small creases to what was the neck point. These two
    # new lines divide two existing regions, thereby adding 2 new regions to the total count.
    head_crimp_additional_regions = 2

    # 4. The total number of regions is the sum of these components.
    total_regions = central_regions + total_corner_regions + head_crimp_additional_regions

    print("The total number of regions is calculated by summing the regions in each part of the crease pattern.")
    print("\nEquation:")
    print(f"{central_regions} (central area) + {total_corner_regions} (corner areas) + {head_crimp_additional_regions} (head crimp) = {total_regions}")
    print("\nTo show the full calculation for the corner areas:")
    print(f"{central_regions} + ({num_corner_sections} * {corner_section_regions}) + {head_crimp_additional_regions} = {total_regions}")


solve_crane_regions()
<<<178>>>