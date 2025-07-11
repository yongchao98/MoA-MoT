def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code based on the specified criteria.
    
    Criteria:
    - Airline: Delta (DL)
    - Stops: 0 (Direct)
    - Fare Class: L, U, T, X, or V (for Platinum Companion Certificate)
    - Format: Standard, minimal length, highest case-insensitive lexicographical sort order.
    """
    
    # Define the components of the extension code
    airline_code = "DL"
    stop_code = "N 0"
    fare_class_code = "f L|U|T|X|V"
    
    # List the components to find the lexicographically highest permutation
    components = [airline_code, stop_code, fare_class_code]
    
    # Sort components in reverse lexicographical order to build the highest-value string
    # Python's default sort is case-sensitive, so we use a key for case-insensitive sort.
    # A reverse sort will place 'N 0' first, then 'f L|U|T|X|V', then 'DL'.
    components.sort(key=str.casefold, reverse=True)
    
    # Join the sorted components with a semicolon to form the final code
    final_code = "; ".join(components)
    
    print(final_code)

generate_ita_matrix_code()