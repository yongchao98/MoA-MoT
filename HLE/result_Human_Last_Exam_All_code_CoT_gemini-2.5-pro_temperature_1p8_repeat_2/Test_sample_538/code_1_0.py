import numpy as np

def calculate_reflections():
    """
    Calculates the number of Bragg reflections for specific plane families
    in a rhombohedrally distorted perovskite (R3m space group) indexed
    on a pseudocubic cell.

    The principle is that the symmetry reduction from cubic to rhombohedral
    causes splitting of certain Bragg peaks. The splitting pattern can be
    predicted by examining the relationship between the plane normals (hkl)
    and the unique rhombohedral axis, which is along the <111> direction
    in the pseudocubic setting.

    Planes that have different angles with respect to the <111> axis will have
    different interplanar spacings (d-spacings) and will thus diffract at
    different angles, resulting in separate observable reflections. We can group
    the planes by calculating the dot product between the plane normal vector
    and the unique axis vector [1, 1, 1]. The number of unique dot product
    magnitudes corresponds to the number of split peaks.
    """
    print("Calculating the number of Bragg reflections for a rhombohedral (R3m) material in a pseudocubic setting.\n")

    # The unique rhombohedral axis in the pseudocubic frame is <111>.
    unique_axis = np.array([1, 1, 1])

    # Define the representative plane normals for each family.
    # These representatives cover all unique orientations within the family.
    plane_families = {
        "{200}": [np.array([2, 0, 0])],
        "{220}": [np.array([2, 2, 0]), np.array([2, -2, 0])],
        "{222}": [np.array([2, 2, 2]), np.array([2, 2, -2])]
    }

    results = {}
    total_reflections = 0

    for family_name, representatives in plane_families.items():
        # Calculate the absolute value of the dot product for each representative normal.
        dot_products = set()
        for normal_vector in representatives:
            # We use the absolute value because the angle itself, not its sign, determines the split.
            # Using int() to avoid float precision issues.
            dot_product_value = abs(np.dot(normal_vector, unique_axis))
            dot_products.add(int(dot_product_value))

        num_reflections = len(dot_products)
        results[family_name] = num_reflections
        total_reflections += num_reflections

        print(f"For the {family_name} family of planes, the number of observed Bragg reflections is: {num_reflections}")

    print("\n--- Summary ---")
    family_counts = list(results.values())
    equation_str = " + ".join(map(str, family_counts))
    print(f"Total Observed Reflections = {equation_str} = {total_reflections}")


if __name__ == "__main__":
    calculate_reflections()