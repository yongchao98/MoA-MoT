def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code based on the specified criteria.
    """
    # Components of the ITA Matrix extension code
    # 1. Fare classes for Platinum Companion Certificate: L, U, T, X, V
    fare_component = "f bc=L|U|T|X|V"
    
    # 2. Airline (Delta marketed and operated), lowercased for optimization
    airline_component = "c:dl"
    
    # 3. Nonstop flights, lowercased for optimization
    stops_component = "nonstop"
    
    # List of components to be sorted
    components = [fare_component, airline_component, stops_component]
    
    # Sort components in descending lexicographical order to maximize the string value
    # str.lower is used for case-insensitive comparison
    components.sort(key=str.lower, reverse=True)
    
    # Join the components with a semicolon to form the final code
    final_code = ";".join(components)
    
    print(final_code)

generate_ita_matrix_code()