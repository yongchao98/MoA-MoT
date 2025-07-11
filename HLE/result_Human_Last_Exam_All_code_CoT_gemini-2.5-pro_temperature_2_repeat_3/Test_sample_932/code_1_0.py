def find_unlikely_tribes():
    """
    Identifies the indices of insect tribes whose immatures are unlikely to be
    collected using a beat-sheet method based on their habitat.
    """
    # Tribe data: (Tribe Name, Immature Habitat)
    # Beat-sheet is effective for collecting insects from external plant foliage.
    tribes = {
        1: ("Apis", "In hive/nest"),
        2: ("Melipotini", "External foliage feeder"),
        3: ("Eupholini", "Internal plant borer"),
        4: ("Acritini", "Dung/carrion/detritus"),
        5: ("Oxyptilini", "External foliage feeder"),
        6: ("Dictyophorini", "External foliage/stem feeder"),
        7: ("Acanthocerini", "Soil/rotten wood")
    }

    unlikely_indices = []
    # Likely habitats are those on the exterior of plants
    likely_habitats = ["External foliage feeder", "External foliage/stem feeder"]

    for index, (name, habitat) in tribes.items():
        if habitat not in likely_habitats:
            unlikely_indices.append(index)

    # Sort the indices in ascending order and format for printing
    unlikely_indices.sort()
    
    # Print the result as a comma-separated string
    result_string = ", ".join(map(str, unlikely_indices))
    print(result_string)

find_unlikely_tribes()