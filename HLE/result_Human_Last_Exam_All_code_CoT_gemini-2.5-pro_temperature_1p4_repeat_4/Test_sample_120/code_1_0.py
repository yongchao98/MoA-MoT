def generate_ita_matrix_code():
    """
    This function generates the ITA Matrix extension code based on the specified criteria.
    """
    # 1. Define the airline code for a direct (single-segment) Delta flight.
    airline_code = "DL"

    # 2. Define the eligible Main Cabin fare classes for a Platinum Delta Companion Certificate.
    fare_classes = ['L', 'U', 'T', 'X', 'V']

    # 3. Sort the fare classes in reverse lexicographical (alphabetical) order
    # to satisfy the "highest value" requirement.
    sorted_fare_classes = sorted(fare_classes, reverse=True)

    # 4. Join the sorted list into a single string separated by '|' and make it lowercase.
    fare_class_string = "|".join(sorted_fare_classes).lower()

    # 5. Assemble the final extension code.
    # The format is ROUTING_CODE;EXTENSION_COMMANDS
    # We remove the space after the semicolon to minimize string length.
    # The extension command for booking code (fare class) is /f bc=
    final_code = f"{airline_code};/f bc={fare_class_string}"

    # 6. Print the final result.
    print(final_code)

generate_ita_matrix_code()