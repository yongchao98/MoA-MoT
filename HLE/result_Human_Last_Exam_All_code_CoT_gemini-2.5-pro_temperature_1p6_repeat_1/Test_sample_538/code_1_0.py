def solve_perovskite_reflections():
    """
    Analyzes the number of Bragg reflections for a rhombohedral (R3m) perovskite
    using pseudocubic indexing.
    """

    # Explanation of the underlying crystallographic principle
    print("This script determines the number of Bragg reflections for a rhombohedral (R3m) material")
    print("when its diffraction pattern is indexed using a pseudocubic cell.")
    print("\nThe transition from a high-symmetry cubic structure to a lower-symmetry rhombohedral")
    print("structure causes certain Bragg peaks to split. A single reflection from a family of planes {hkl}")
    print("in the cubic system splits into multiple reflections if the planes in that family are no longer")
    print("symmetrically equivalent in the lower-symmetry rhombohedral system.")
    print("-" * 75)

    # Analysis for the {200} family of planes
    num_reflections_200 = 1
    print("Analysis for the {200} family of planes:")
    print("In the parent cubic structure, the {200} family includes planes like (200), (020), and (002).")
    print("In the rhombohedral system (Laue group -3m), the 3-fold rotational symmetry that")
    print("relates these planes is preserved. Therefore, they remain symmetrically equivalent, have the")
    print("same d-spacing, and contribute to a single Bragg reflection.")
    print("\nResult for {200}:")
    print(f"Number of observed Bragg reflections for {{200}} = {num_reflections_200}")
    print("-" * 75)

    # Analysis for the {220} family of planes
    num_reflections_220 = 2
    print("Analysis for the {220} family of planes:")
    print("In the cubic system, the {220} family includes planes like (220), (2-20), (202), etc.")
    print("When symmetry is lowered to rhombohedral, these planes are no longer all equivalent.")
    print("They split into two distinct groups based on their orientation:")
    print("  1. Planes of the type (hh0), such as (220), (202), and (022).")
    print("  2. Planes of the type (h-h0), such as (2-20), (20-2), and (0-22).")
    print("These two groups have different d-spacings and thus produce two separate Bragg reflections.")
    print("\nResult for {220}:")
    print(f"Number of observed Bragg reflections for {{220}} = {num_reflections_220}")
    print("-" * 75)

    # Analysis for the {222} family of planes
    num_reflections_222 = 1
    print("Analysis for the {222} family of planes:")
    print("The {222} planes are perpendicular to the <111> body diagonals of the cubic cell.")
    print("Since the rhombohedral distortion occurs along one of these <111> directions,")
    print("this particular set of planes retains its high symmetry.")
    print("All planes in the {222} family remain equivalent and produce a single, unsplit reflection.")
    print("\nResult for {222}:")
    print(f"Number of observed Bragg reflections for {{222}} = {num_reflections_222}")
    print("-" * 75)

# Execute the analysis
solve_perovskite_reflections()
<<<1, 2, 1>>>