def generate_ita_matrix_code():
    """
    This function constructs and prints the optimized ITA Matrix extension code.
    """

    # --- Component Breakdown ---
    # 'o:' specifies the outbound portion of the trip.
    # 'dl' is the airline code for Delta.
    # We use lowercase as permitted by the prompt.
    carrier_and_leg = "o:dl"

    # '/-s' is the shortest qualifier for a non-stop flight.
    # It is a routing rule, so it comes before the semicolon.
    routing_rules = "/-s"

    # ';' separates routing rules from fare rules.
    separator = ";"

    # 'f' is the shorthand for fare-related rules.
    fare_rule_qualifier = "f"

    # Delta Platinum Companion Certificates are valid for fare classes L, U, T, X, V.
    # To meet the optimization requirement of highest lexicographical value,
    # we sort these characters in reverse alphabetical order.
    # We use lowercase as permitted.
    fare_codes = "x|v|u|t|l"

    # --- Final Assembly ---
    # We combine the components, removing spaces to minimize string length.
    # The final code is structured as: [leg]:[carrier][routing];[fare_qualifier] [fare_codes]
    final_code = f"{carrier_and_leg}{routing_rules}{separator}{fare_rule_qualifier} {fare_codes}"
    
    # As ITA is space-tolerant and removing the space after 'f' creates the highest
    # lexicographical value for minimal length strings, we will use that version.
    # Let's compare "f x..." and "fx..."
    # Lexicographically, ' ' (space) comes before 'x', so 'fx...' is greater than 'f x...'.
    # Therefore, we remove the space.
    # Minimal length: o:dl/-s;fl|u|t|v|x -> has length 18
    # Sorting again: o:dl/-s;fx|v|u|t|l is the best candidate.

    # Re-evaluating minimal length and sorting.
    # ITA Matrix syntax examples often show `f L|U|T..`. A space is not required.
    # So `o:dl/-s;fx|v|u|t|l` is a valid, 18-character string.
    # Let's consider `o:dl/-s;f l|u..`. The space makes it longer.
    # My previous plan `o:dl/-s;f x|v|u|t|l` (with space) is 19 characters.
    # The spaceless version `o:dl/-s;fx|v|u|t|l` is 18 characters and thus minimal length.
    # This string also has the highest lexicographical value for its length.
    final_code_minimal = f"{carrier_and_leg}{routing_rules}{separator}{fare_rule_qualifier}{fare_codes}"


    print("The final ITA Matrix extension code is:")
    print(final_code_minimal)


generate_ita_matrix_code()