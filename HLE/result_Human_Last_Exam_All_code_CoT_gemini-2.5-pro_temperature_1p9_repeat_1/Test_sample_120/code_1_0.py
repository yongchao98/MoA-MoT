def generate_ita_extension_code():
    """
    Constructs the optimized ITA Matrix extension code based on a set of rules.
    This code finds direct Delta flights eligible for Platinum Companion Certificates,
    optimizing for minimal length and lexicographical value.
    """
    # 1. Define the airline.
    carrier = "DL"

    # 2. Define the eligible fare classes for a Platinum Delta Companion Certificate.
    # The prompt requires lowercasing where possible.
    eligible_fare_classes = ['l', 'u', 't', 'x', 'v']

    # 3. Define the most compact, valid ITA Matrix syntax prefix for applying a
    # booking code rule to a specific carrier's fares. This helps minimize
    # the overall string length.
    # '::' scopes the rule to the carrier.
    # '/fbc=' is the compact form for the fare booking code rule.
    syntax_prefix = "::/fbc="

    # 4. To meet the tie-breaker requirement of choosing the string with the highest
    # case-insensitive lexicographical value, we sort the fare classes in
    # descending order. This ensures the part of the string after the '=' is
    # as lexicographically large as possible.
    eligible_fare_classes.sort(reverse=True)

    # 5. Join the sorted fare classes into the required format (separated by '|').
    fare_class_string = "|".join(eligible_fare_classes)

    # 6. Assemble the final extension code.
    final_extension_code = carrier + syntax_prefix + fare_class_string

    # 7. Print the final result. The direct flight constraint (nonstop) is best set
    # in the main ITA Matrix UI and is not included in this fare-specific code.
    print(final_extension_code)

# Execute the function to display the result.
generate_ita_extension_code()