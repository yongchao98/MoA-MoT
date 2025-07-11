def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code for finding Delta direct flights
    eligible for Platinum Companion Certificates, following specific formatting rules.
    """

    # 1. Define the components of the extension code.
    # Specifier for a direct (0-stop) outbound flight on Delta.
    flight_specifier = "O:DL"

    # Extension code command to filter by fare (booking class).
    # Lowercased as requested.
    fare_qualifier = "/f"

    # Eligible fare classes for the Platinum Delta Companion Certificate.
    fare_classes = ['L', 'U', 'T', 'X', 'V']

    # 2. Sort fare classes in descending order to maximize lexicographical value.
    # This addresses the tie-breaker rule.
    sorted_fare_classes = sorted(fare_classes, reverse=True)

    # 3. Join the sorted fare classes with the '|' (OR) operator.
    fare_class_string = "|".join(sorted_fare_classes)

    # 4. Assemble the final extension code string.
    final_code = f"{flight_specifier} {fare_qualifier} {fare_class_string}"

    # Print each part of the "equation" as requested.
    print("Constructing the ITA Matrix Extension Code:")
    print(f"Part 1 (Outbound Direct Flight on Delta): {flight_specifier}")
    print(f"Part 2 (Fare Class Filter command): {fare_qualifier}")
    print(f"Part 3 (Eligible Fare Classes, sorted for max lexicographical value): {fare_class_string}")
    print("\nFinal Resulting Code:")
    print(final_code)


if __name__ == "__main__":
    generate_ita_matrix_code()
