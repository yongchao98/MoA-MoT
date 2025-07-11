def generate_ita_extension_code():
    """
    Constructs the ITA Matrix extension code for finding direct Delta flights
    eligible for Platinum Companion Certificates, following specific optimization rules.
    """

    # Rule: Direct (non-stop) flights on Delta (dl).
    # The 'O:' command is the shortest syntax for non-stop.
    airline_and_stops_command = "O:dl"

    # Rule: Platinum Companion Certificate eligible fare classes are L, U, T, X, and V.
    # The 'C' command specifies fare classes. Lowercasing and removing spaces minimizes length.
    fare_class_command = "Cl,u,t,x,v"

    # The two primary components of the code are combined with '::'.
    # To satisfy the "lexicographically highest" rule among strings of minimal length,
    # we must compare the two possible orderings.
    option1 = f"{airline_and_stops_command}::{fare_class_command}"
    option2 = f"{fare_class_command}::{airline_and_stops_command}"

    # 'O' comes after 'C' in the alphabet, so option1 is lexicographically higher.
    # We select the lexicographically highest valid code of minimal length.
    final_code = max(option1, option2)

    print(final_code)

generate_ita_extension_code()