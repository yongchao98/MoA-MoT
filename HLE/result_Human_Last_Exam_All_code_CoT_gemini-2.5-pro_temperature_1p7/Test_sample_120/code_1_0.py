def solve_ita_matrix_code():
    """
    This script determines the ITA Matrix extension code for finding direct Delta
    flights eligible for Platinum Companion Certificates, applying specified
    optimization rules.
    """

    # 1. Define the routing code for a direct Delta flight.
    # The shortest syntax is the airline's two-letter code.
    # Per the instructions, we use lowercase.
    routing_code = "dl"

    # 2. Define the eligible fare classes for the Platinum Companion Certificate.
    # These must be uppercase as they are standard fare codes.
    fare_classes = ["L", "U", "T", "X", "V"]

    # 3. Sort the fare classes in descending alphabetical order.
    # This is to satisfy the "highest value in case-insensitive lexicographic sorting"
    # rule for strings of minimal length.
    fare_classes.sort(reverse=True)
    # The sorted list will be ['X', 'V', 'U', 'T', 'L']

    # 4. Join the fare classes into a single string using the "|" OR operator.
    fare_class_string = "|".join(fare_classes)

    # 5. Assemble the final extension code string.
    # The format is <routing> ; f bc=<fares>
    # 'f' is the fare qualifier command, and 'bc' specifies the booking code (fare class).
    final_code = f"{routing_code} ; f bc={fare_class_string}"

    # 6. Print the final result.
    print(final_code)

solve_ita_matrix_code()