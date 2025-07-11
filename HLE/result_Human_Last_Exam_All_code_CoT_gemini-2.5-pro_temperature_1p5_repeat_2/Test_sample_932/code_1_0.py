def solve_insect_collection():
    """
    Determines which insect tribes are unlikely to be collected with a beat-sheet method
    by analyzing the typical habitat of their immature stages.
    """

    # A beat-sheet is used to collect insects living externally on plant foliage.
    # We identify tribes whose immatures live in concealed or different habitats.
    
    # Format: {index: (Tribe Name, Immature Habitat, Collectible?)}
    tribes_data = {
        1: ("Apis", "Inside a hive/nest", False),
        2: ("Melipotini", "Externally on foliage (caterpillars)", True),
        3: ("Eupholini", "Inside plant tissue (wood/stem boring grubs)", False),
        4: ("Acritini", "In dung, carrion, or under bark (predatory larvae)", False),
        5: ("Oxyptilini", "Externally on foliage (caterpillars)", True),
        6: ("Dictyophorini", "Externally on foliage (nymphs)", True),
        7: ("Acanthocerini", "In soil or rotting wood (grubs)", False)
    }

    unlikely_indices = []
    for index, data in tribes_data.items():
        is_collectible = data[2]
        if not is_collectible:
            unlikely_indices.append(index)
    
    # Sort the indices in ascending order
    unlikely_indices.sort()
    
    # The final answer is the indices of the unlikely tribes, sorted and comma-separated.
    result_string = ",".join(map(str, unlikely_indices))
    print(result_string)

solve_insect_collection()