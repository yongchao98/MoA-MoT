def generate_ita_matrix_code():
    """
    This function generates an ITA Matrix extension code based on a specific set of rules.
    """
    # 1. Define the individual components of the extension code.
    # - Direct flights (0 connections). Command is 'maxconnect'.
    # - Delta as the marketing carrier. Command is 'c'.
    # - Eligible Platinum Companion Certificate fare classes (L, U, T, X, V).
    #   Command is 'f' for fare rules, 'bc' for booking class.
    #
    # Components are lowercased where possible per the prompt.
    components = [
        "maxconnect 0",
        "f bc=l|u|t|x|v",
        "c:DL"
    ]

    # 2. To find the string with the highest case-insensitive lexicographical value,
    #    sort the components in reverse alphabetical order.
    components.sort(key=str.lower, reverse=True)

    # 3. Join the sorted components with a semicolon.
    #    No spaces are used to ensure minimal string length.
    final_code = ";".join(components)

    # 4. Print the final result.
    print(final_code)

generate_ita_matrix_code()