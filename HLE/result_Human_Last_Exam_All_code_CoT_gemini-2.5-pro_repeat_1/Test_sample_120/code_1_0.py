def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code based on the specified rules.
    """
    # 1. Routing code for a direct (non-stop) flight on Delta.
    routing_code = "DL"

    # 2. Fare classes for the Platinum Delta Companion Certificate.
    fare_classes = ["L", "U", "T", "X", "V"]
    
    # 3. Sort the fare classes alphabetically for a deterministic result.
    sorted_fare_classes = sorted(fare_classes)
    
    # 4. Join the sorted classes with the '|' (OR) operator.
    fare_class_string = "|".join(sorted_fare_classes)

    # 5. Fare class command, lowercased as per instructions.
    fare_command = "f"

    # 6. Separator. We choose ';' over '/' because it has a higher
    # lexicographical value, as required by the tie-breaking rule.
    separator = ";"

    # 7. Assemble the final code. A space is required between the
    # fare command 'f' and the list of fare classes.
    final_code = f"{routing_code}{separator}{fare_command} {fare_class_string}"

    # Print the "equation" showing how the code was built.
    # The prompt asks to output each part of the final equation.
    print("Building the extension code:")
    print(f"Routing Code (Direct Delta): {routing_code}")
    print(f"Separator (Lexicographically Highest): '{separator}'")
    print(f"Fare Command (Lowercase): {fare_command}")
    print(f"Eligible Fare Classes (Sorted): {fare_class_string}")
    print("-" * 20)
    print("Final Extension Code:")
    print(final_code)
    
    # Hidden final answer block for the system.
    print(f"\n<<<__{final_code}__>>>")

generate_ita_matrix_code()