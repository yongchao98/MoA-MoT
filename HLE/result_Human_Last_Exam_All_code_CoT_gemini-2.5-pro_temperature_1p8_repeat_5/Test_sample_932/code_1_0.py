def find_unlikely_tribes():
    """
    Identifies insect tribes whose immatures are unlikely to be collected
    using a beat-sheet method based on their life habits.
    """

    # Tribe names corresponding to the indices
    tribes = {
        1: "Apis",
        2: "Melipotini",
        3: "Eupholini",
        4: "Acritini",
        5: "Oxyptilini",
        6: "Dictyophorini",
        7: "Acanthocerini"
    }

    # Life habits of the immatures and suitability for beat-sheet collection
    # True = Unlikely to be collected, False = Likely to be collected
    unlikely_habits = {
        1: "Immatures (larvae) live protected inside a hive, not on plants.",
        2: "Immatures (caterpillars) are external foliage feeders.",
        3: "Immatures (larvae) are typically internal borers or live in soil.",
        4: "Immatures (larvae) are predators in dung, carrion, or leaf litter, not on plants.",
        5: "Immatures (caterpillars) are external feeders on plant parts.",
        6: "Immatures (nymphs) live and feed externally on plants.",
        7: "Immatures (larvae/grubs) live in soil or rotting wood."
    }
    
    unlikely_indices = []
    
    print("Analyzing tribes for beat-sheet collection suitability of immatures:\n")

    for index in sorted(tribes.keys()):
        tribe_name = tribes[index]
        reason = unlikely_habits[index]
        # Check if the reason indicates an uncollectible habit
        if "not on plants" in reason or "internal" in reason or "in soil" in reason or "in hive" in reason or "in dung" in reason or "in rotting wood" in reason:
            print(f"Index {index} ({tribe_name}): Unlikely. Reason: {reason}")
            unlikely_indices.append(index)
        else:
            print(f"Index {index} ({tribe_name}): Likely. Reason: {reason}")
            
    # Sort the indices in ascending order
    unlikely_indices.sort()
    
    # Format the final answer as a comma-separated string
    final_answer_string = ", ".join(map(str, unlikely_indices))
    
    print("\nThe indices of the tribes whose immatures are unlikely to be collected, in ascending order, are:")
    print(final_answer_string)

find_unlikely_tribes()