def generate_ita_matrix_code():
    """
    Constructs the optimized ITA Matrix extension code based on the specified rules.
    """
    # 1. Define the components based on the search requirements.
    # Carrier and direct flight requirement. 'dl' is the shortest routing code for a
    # single direct flight on Delta (marketing carrier).
    routing_code = "dl"

    # Fare rule qualifier for booking codes.
    fare_rule_prefix = "f bc="

    # Eligible booking codes for a Platinum Delta Companion Certificate.
    booking_codes = ['L', 'U', 'T', 'X', 'V']

    # 2. Apply optimization rules.
    # Rule: Use lowercase where possible.
    # We will use lowercase for all parts of the string.
    booking_codes_lower = [code.lower() for code in booking_codes]

    # Rule: Maximize case-insensitive lexicographical value.
    # To do this, sort the booking codes in reverse alphabetical order.
    booking_codes_sorted = sorted(booking_codes_lower, reverse=True)

    # Join the sorted codes with the '|' (OR) operator.
    fare_codes_string = "|".join(booking_codes_sorted)

    # Rule: Minimize length and maximize lexicographical value for the separator.
    # A space ' ' and a semicolon ';' are both valid separators.
    # 'dl;f...' and 'dl f...' have the same minimal length.
    # The ASCII value of ';' is higher than ' ', so we choose the semicolon
    # to maximize the lexicographical value of the overall string.
    separator = ";"

    # 3. Assemble the final extension code string.
    final_code = routing_code + separator + fare_rule_prefix + fare_codes_string

    print("ITA Matrix Outbound Extension Code:")
    print("This code specifies a direct flight on Delta ('dl') with fare booking codes")
    print(f"limited to {', '.join(booking_codes)} for Companion Certificate eligibility.")
    print("The components are optimized for minimal length and highest lexicographical value.")
    print("\nFinal Code:")
    print(final_code)

if __name__ == "__main__":
    generate_ita_matrix_code()