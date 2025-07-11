def generate_ita_matrix_code():
    """
    This function generates the ITA Matrix extension code for finding direct
    Delta flights eligible for Platinum Companion Certificates, optimized for
    length and lexicographical order.
    """
    # Components of the extension code
    # 1. Direct flight (0 stops)
    direct_flight_rule = "stops 0"
    
    # 2. Fare class rule for Platinum Companion Certificate (L, U, T, X, V classes)
    fare_class_rule = "f bc=L|U|T|X|V"
    
    # 3. Airline rule for Delta
    airline_rule = "DL"
    
    # List of components to sort
    components = [direct_flight_rule, fare_class_rule, airline_rule]
    
    # Sort components in descending lexicographical order to get the highest value string
    components.sort(reverse=True)
    
    # Join the components with the shortest separator ';'
    final_code = "; ".join(components)
    
    print(final_code)

generate_ita_matrix_code()