def generate_ita_matrix_code():
    """
    This function generates the ITA Matrix extension code for finding
    direct Delta flights eligible for Platinum Companion Certificates.
    """

    # 1. Define the eligible fare classes for Platinum Companion Certificates.
    fare_classes = ['L', 'U', 'T', 'X', 'V']

    # 2. Per the prompt, use lowercase where possible.
    lower_fare_classes = [c.lower() for c in fare_classes]

    # 3. Sort fare classes in reverse lexicographical order to maximize the
    #    string's sort value, as required by the tie-breaking rule.
    sorted_classes = sorted(lower_fare_classes, reverse=True)

    # 4. Construct the booking code (bc) filter string.
    #    The format is bc=CLASS1|bc=CLASS2|...
    booking_code_filter_parts = [f"bc={code}" for code in sorted_classes]
    booking_code_filter = "|".join(booking_code_filter_parts)

    # 5. Define other components of the extension code.
    #    - 'f' specifies a fare rule, which is lexicographically higher than 'DL'.
    #    - 'DL' is the carrier code for Delta.
    #    - '/nonstop' is the shortest modifier for direct flights.
    fare_rule_prefix = "f"
    carrier_and_stop_rule = "DL /nonstop"

    # 6. Assemble the final code, combining the fare and routing rules.
    #    This structure was chosen based on the length and lexicographical sort rules.
    final_extension_code = f"{fare_rule_prefix} {booking_code_filter}; {carrier_and_stop_rule}"

    print(final_extension_code)

if __name__ == "__main__":
    generate_ita_matrix_code()