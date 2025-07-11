def solve_ita_matrix_code():
    """
    This function determines the optimal ITA Matrix extension code based on a set of rules.
    """
    
    # 1. Define the components of the extension code.
    airline_and_stops_code = "DL"  # Specifies a single, direct flight on Delta.
    
    # Fare classes eligible for Platinum Delta Companion Certificates.
    eligible_fare_classes = ['L', 'U', 'T', 'X', 'V']
    
    # 2. Apply optimization rules.
    
    # Rule: Maximize case-insensitive lexicographical value.
    # To do this, we sort the fare classes in reverse alphabetical order.
    eligible_fare_classes.sort(reverse=True)
    
    # Join the sorted fare classes with the '|' (OR) operator.
    fare_class_string = "|".join(eligible_fare_classes)
    
    # Rule: Minimize overall string length.
    # We choose the shortest valid syntax.
    # 'fbc=' is shorter than 'f bc='.
    # Between the routing/fare separators ';' and '/', both are one character,
    # so length is the same.
    #
    # Rule: Among minimal length strings, choose the one with the highest
    # lexicographical value.
    # Comparing ';' (ASCII 59) and '/' (ASCII 47), ';' is higher.
    separator = ";"
    
    # The fare rule portion of the code.
    # 'f' and 'bc' are lowercased as requested.
    fare_rule_code = f"fbc={fare_class_string}"
    
    # 3. Assemble the final extension code.
    final_extension_code = f"{airline_and_stops_code}{separator}{fare_rule_code}"
    
    print(final_extension_code)

solve_ita_matrix_code()