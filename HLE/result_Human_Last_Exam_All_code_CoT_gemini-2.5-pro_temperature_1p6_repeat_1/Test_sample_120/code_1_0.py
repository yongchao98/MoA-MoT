def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code based on the specified rules.
    """
    # 1. Define the components
    # Airline code for Delta. Lowercased as per instruction.
    airline_code = "dl"

    # Eligible fare classes for Platinum Companion Certificate. Lowercased as per instruction.
    fare_classes = ['l', 'u', 't', 'x', 'v']

    # 2. Sort the fare classes for the highest lexicographical value
    # Sorting in reverse (descending) order.
    fare_classes.sort(reverse=True)

    # 3. Construct the fare class specification string
    # e.g., "bc=x|bc=v|bc=u|bc=t|bc=l"
    fare_class_parts = []
    for fc in fare_classes:
        fare_class_parts.append(f"bc={fc}")
    fare_string = "|".join(fare_class_parts)

    # 4. Assemble the final extension code
    # The format is: [airline] /f [fare_string]
    # Example: "dl /f bc=x|bc=v|bc=u|bc=t|bc=l"
    final_code = f"{airline_code} /f {fare_string}"

    # Print the final result
    print(final_code)

generate_ita_matrix_code()