import sys

def solve():
    """
    This script determines the correct ITA Matrix extension code based on a specific set of rules.
    It constructs the code step-by-step and then prints the final, optimized result.
    """

    # 1. Airline: Delta (DL). Lowercased for lexicographical sorting tie-breaker ('d' > 'D').
    airline_code = "dl"

    # 2. Flight Type: Direct (no layover) means non-stop.
    # The shortest qualifier is '/stops 0'.
    non_stop_qualifier = "/stops 0"

    # 3. Fare Classes: Platinum Companion Certificate is valid for L, U, T, X, V.
    # To maximize lexicographical value, sort these in reverse alphabetical order.
    fare_codes = sorted(['L', 'U', 'T', 'X', 'V'], reverse=True)
    fare_codes_str = ",".join(fare_codes)
    fare_class_qualifier = f"/f bc={fare_codes_str}"

    # 4. Assemble and Optimize for Lexicographical Value:
    # Qualifiers should be sorted in reverse alphabetical order to maximize the string's value.
    # '/stops' comes after '/f', so it goes first.
    qualifiers = [non_stop_qualifier, fare_class_qualifier]
    
    # Sort qualifiers based on their first letter after the slash, in reverse order
    # Example: '/stops' (s) vs '/f' (f) -> '/stops' comes first
    qualifiers.sort(key=lambda x: x[1], reverse=True)

    # 5. Combine all parts into the final extension code string.
    # The structure is: [airline] [qualifier1] [qualifier2] ...
    final_extension_code = f"{airline_code} {' '.join(qualifiers)}"

    # Print the final result
    print(final_extension_code)

solve()