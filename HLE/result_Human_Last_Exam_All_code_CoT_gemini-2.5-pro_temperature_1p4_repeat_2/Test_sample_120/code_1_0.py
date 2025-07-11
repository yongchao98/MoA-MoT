import sys

def generate_ita_matrix_code():
    """
    Generates the ITA Matrix extension code for finding direct Delta flights
    eligible for a Platinum Companion Certificate, optimized for length and
    lexicographical order.
    """
    # 1. Define the components of the search query.
    # Airline specification for Delta.
    carrier_specifier = "C DL"
    
    # Eligible fare classes for the Platinum Companion Certificate.
    fare_classes = ['L', 'U', 'T', 'X', 'V']
    
    # 2. Prepare commands based on optimization rules.
    # The command for direct flights.
    direct_flight_command = "nonstop"
    
    # Format the fare class command: lowercase and sorted alphabetically.
    # This ensures the output is deterministic and helps with the final sorting rule.
    sorted_fare_classes = sorted([f.lower() for f in fare_classes])
    fare_class_command = f"f bc={'|'.join(sorted_fare_classes)}"
    
    # 3. Assemble the extension code parts.
    # The two main commands to be chained together.
    commands = [direct_flight_command, fare_class_command]
    
    # 4. Apply the tie-breaker: choose the highest value in case-insensitive
    # lexicographic sorting.
    # Python's sort (or sorted) on strings is lexicographical. We want the
    # highest value, so we sort in reverse order and take the first element.
    commands.sort(key=str.lower, reverse=True)
    
    # 5. Construct the final extension code string.
    # The format is C <carrier>;/<command1>;<command2>...
    final_code = f"{carrier_specifier};/{';'.join(commands)}"
    
    print(final_code)

generate_ita_matrix_code()