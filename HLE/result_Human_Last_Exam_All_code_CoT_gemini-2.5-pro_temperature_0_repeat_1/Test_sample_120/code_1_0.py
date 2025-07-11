def generate_ita_matrix_code():
    """
    This script determines the ITA Matrix extension code for finding direct Delta flights
    eligible for Platinum Companion Certificates, optimized for length and lexicographical value.
    """

    # 1. Define the components for the search query.
    # Component for direct flights (0 stops).
    direct_flights_code = "maxstops 0"
    
    # Component for Platinum Companion Certificate fare classes (L, U, T, X, V).
    fare_class_code = "f bc=l|u|t|x|v"
    
    # Component for specifying Delta as the marketing carrier.
    airline_code = "c:dl"

    # 2. Explain each component of the final code.
    print("The final extension code is constructed from the following components:")
    
    print(f"\n1. To limit the search to direct flights, we specify a maximum of 0 stops.")
    print(f"   - Code: {direct_flights_code}")
    
    print(f"\n2. To filter for fares eligible for the Platinum Companion Certificate, we specify the allowed booking codes (L, U, T, X, V).")
    print(f"   - Code: {fare_class_code}")

    print(f"\n3. To ensure the flight is a Delta flight, we specify the marketing carrier.")
    print(f"   - Code: {airline_code}")

    # 3. Combine the components based on the optimization rules.
    # To achieve the highest case-insensitive lexicographical value, the components are
    # ordered as follows: 'maxstops...' > 'f bc...' > 'c:dl'.
    components_ordered = [direct_flights_code, fare_class_code, airline_code]
    final_code = ";".join(components_ordered)

    print("\nThese components are ordered to create the valid string of minimal length with the highest lexicographical value.")
    print("\nFinal Outbound Extension Code:")
    print(final_code)
    
    # The final answer is also provided in the requested format.
    print(f"\n<<<{final_code}>>>")

generate_ita_matrix_code()