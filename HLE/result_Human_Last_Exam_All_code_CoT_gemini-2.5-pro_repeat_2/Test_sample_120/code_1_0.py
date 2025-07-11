def generate_ita_matrix_code():
    """
    Constructs the ITA Matrix extension code based on a set of rules for a
    specific Delta flight search.
    """
    # 1. Define the components of the extension code based on the requirements.
    #    - Airline: Delta (DL)
    #    - Direct Flight: 0 connections (X 0)
    #    - Fare Classes: Platinum Companion Cert (L, U, T, X, V)
    #    Commands are lowercased as per the prompt's formatting rule.
    airline_component = "c DL"
    stops_component = "x 0"
    fares_component = "f L|U|T|X|V"

    # 2. Put the components into a list for sorting.
    components = [airline_component, stops_component, fares_component]

    # 3. Sort the components in reverse (descending) lexicographical order
    #    to satisfy the "highest value" requirement. String comparison handles this
    #    correctly ('x' > 'f' > 'c').
    components.sort(reverse=True)

    # 4. Join the sorted components with a semicolon to form the final code.
    final_extension_code = ";".join(components)

    # 5. Print the final result.
    print(final_extension_code)

generate_ita_matrix_code()