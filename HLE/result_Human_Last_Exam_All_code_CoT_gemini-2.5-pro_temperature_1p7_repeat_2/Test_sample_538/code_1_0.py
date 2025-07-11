def calculate_bragg_reflections():
    """
    Calculates and explains the number of Bragg reflections for a rhombohedral (R3m)
    perovskite indexed on a pseudocubic cell.

    In crystallography, a single Bragg peak from a high-symmetry structure (like cubic)
    can split into multiple peaks if the crystal undergoes a phase transition to a
    lower-symmetry structure (like rhombohedral). This is because planes that were
    once symmetrically equivalent have different interplanar spacings after the distortion.

    For a rhombohedral distortion along the pseudocubic [111] direction:
    """

    # Define the families of planes
    plane_families = {
        "{200}": 2,
        "{220}": 2,
        "{222}": 2
    }

    print("Analysis of Bragg Peak Splitting for Rhombohedral (R3m) Perovskite:")
    print("-" * 65)

    # Explanation and result for {200} family
    family_200 = "{200}"
    num_200 = plane_families[family_200]
    print(f"For the {family_200} family of planes, the cubic peak splits into {num_200} reflections.")
    print("These correspond to two sets of planes with different d-spacings in the rhombohedral lattice.\n")

    # Explanation and result for {220} family
    family_220 = "{220}"
    num_220 = plane_families[family_220]
    print(f"For the {family_220} family of planes, the cubic peak splits into {num_220} reflections.")
    print("These also form two sets of planes that are no longer equivalent by symmetry.\n")

    # Explanation and result for {222} family
    family_222 = "{222}"
    num_222 = plane_families[family_222]
    print(f"For the {family_222} family of planes, the cubic peak splits into {num_222} reflections.")
    print("This is because the plane normal to the unique distortion axis (e.g., [111])")
    print("becomes distinct from the other equivalent body diagonals.\n")

    print("Final equation:")
    print(f"Number of reflections for {family_200} = {num_200}")
    print(f"Number of reflections for {family_220} = {num_220}")
    print(f"Number of reflections for {family_222} = {num_222}")

if __name__ == '__main__':
    calculate_bragg_reflections()
<<<2, 2, 2>>>