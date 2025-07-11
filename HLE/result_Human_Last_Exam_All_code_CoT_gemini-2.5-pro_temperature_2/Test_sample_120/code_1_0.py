import sys

def solve():
    """
    This script determines the correct ITA Matrix extension code for finding direct
    Delta flights eligible for a Platinum Companion Certificate, optimized for
    minimum length and maximum lexicographical value.
    """

    # 1. Define the fare classes for a Platinum Companion Certificate.
    # As per the prompt, they must be lowercased.
    fare_classes = ['l', 'u', 't', 'x', 'v']

    # 2. To maximize lexicographical value, sort the classes in reverse order.
    # This ensures that when joined into a string, it's as "high" as possible
    # alphabetically (e.g., "x" comes before "v").
    fare_classes.sort(reverse=True)
    fare_class_string = ",".join(fare_classes)

    # 3. Define the other extension code components.
    # 'DL' specifies the airline (Delta).
    # 'C0' is the shortest way to specify 0 connections (a direct flight).
    # The '0' is the number mentioned in the prompt.
    airline_part = f"DL /fclass {fare_class_string}"
    connection_part = "C0"

    # 4. We have two parts to combine: airline/fare and connection.
    # To maximize the final string's lexicographical value, they must be
    # ordered from highest to lowest. Comparing them, "DL..." is greater than "C0".
    code_parts = [airline_part, connection_part]

    # In this specific case, we know 'DL...' > 'C0' but sorting ensures correctness.
    code_parts.sort(key=lambda x: x.lower(), reverse=True)

    # 5. Join the optimally ordered parts with a semicolon.
    final_extension_code = ";".join(code_parts)

    print(final_extension_code)

solve()