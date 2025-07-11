def generate_ita_matrix_code():
    """
    Constructs and prints the ITA Matrix extension code based on the specified criteria.
    """
    # Base airline code for Delta
    airline = "DL"

    # Rule for direct (non-stop) flights
    # 'f' is a flight-specific rule. 'bc=0' specifies 0 connections.
    direct_flight_rule = "f bc=0"

    # Rule for Platinum Companion Certificate eligible fare classes
    # '/f' is a fare-specific rule. 'bc=' specifies the booking class.
    # The fare classes (L, U, T, X, V) are sorted in reverse alphabetical
    # order (X, V, U, T, L) to maximize the lexicographical value of the final string.
    fare_class_rule = "/f bc=x|v|u|t|l"

    # To maximize the case-insensitive lexicographical value of the string,
    # the rules are ordered based on their first character ('f' > '/').
    # The separator '::' is used between rules.
    final_code = f"{airline}::{direct_flight_rule}::{fare_class_rule}"

    print(final_code)

generate_ita_matrix_code()