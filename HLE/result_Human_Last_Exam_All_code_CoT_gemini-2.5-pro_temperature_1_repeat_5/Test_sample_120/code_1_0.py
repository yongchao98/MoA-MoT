def generate_ita_matrix_code():
    """
    Constructs and prints the optimized ITA Matrix extension code based on the given constraints.
    """
    # Components of the extension code
    # 1. Airline: Delta (DL). Must be uppercase.
    airline_code = "DL"

    # 2. Non-stop flight. The shortest command is 'N'.
    non_stop_code = "N"

    # 3. Fare classes for Platinum Companion Certificate: L, U, T, X, V.
    # To maximize lexicographical value, we sort them in descending order.
    # We also use lowercase as requested.
    fare_classes = sorted(['l', 'u', 't', 'x', 'v'], reverse=True)
    fare_class_string = "|".join(fare_classes)
    fare_class_code = f"f bc={fare_class_string}"

    # 4. Outbound/Inbound separator for outbound-only rules.
    outbound_separator = "::"

    # Assemble the components. To maximize lexicographical value, we sort the
    # individual rule components ('N', 'f bc=...', 'DL') in descending
    # case-insensitive alphabetical order.
    # N > f > d
    components = [non_stop_code, fare_class_code, airline_code]

    # Join with the ITA Matrix rule separator ';'
    final_code = "; ".join(components) + f" {outbound_separator}"

    print(final_code)

generate_ita_matrix_code()
#<<<N; f bc=x|v|u|t|l; DL ::>>>