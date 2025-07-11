def generate_ita_matrix_code():
    """
    This function determines the ITA Matrix extension code for finding flights
    eligible for a Platinum Delta Companion Certificate, following all formatting rules.
    """

    # --- Component Definitions ---

    # 1. Fare Classes: Platinum Companion Certificates are valid for L, U, T, X, V.
    # We use the concise 'f' qualifier and lowercase the fare codes.
    fare_component = "f l|u|t|x|v"

    # 2. Operating Carrier: The certificate requires flights to be Delta-operated.
    # To minimize the extension code's length, we assume the user specifies the
    # marketing carrier (DL) and non-stop preference in the main routing fields.
    # This leaves the operator requirement for the extension code.
    operator_component = "operator dl"

    # --- Assembly and Formatting ---

    # List of components to be included in the extension code.
    components = [fare_component, operator_component]

    # Sort the components in reverse (descending) lexicographical order
    # to produce the final string with the highest value.
    components.sort(key=str.casefold, reverse=True)

    # Join the components with the standard ITA Matrix separator.
    final_code = " ; ".join(components)

    print(final_code)

generate_ita_matrix_code()