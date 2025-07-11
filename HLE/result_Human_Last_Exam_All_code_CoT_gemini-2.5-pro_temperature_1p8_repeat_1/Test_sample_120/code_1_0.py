import itertools

def solve_ita_matrix_code():
    """
    This script determines the optimal ITA Matrix extension code based on a specific set of rules.
    """
    print("Finding the ITA Matrix extension code for direct Delta flights eligible for Platinum Companion Certificates.")
    print("The final code must have the minimum possible length and, among those, the highest case-insensitive lexicographical value.\n")

    # 1. Define the components of the extension code using their shortest form.

    # To specify the carrier is Delta (DL), the shortest command is "f dl".
    # This is shorter than "/airline dl" or "carrier dl".
    airline_component = "f dl"

    # To specify a direct (non-stop) flight, the command is "maxconnect 0".
    # The number '0' here explicitly sets the maximum number of layovers to zero.
    stops_component = "maxconnect 0"

    # To filter for fares eligible for a Platinum Companion Certificate (main cabin),
    # we need to match fare classes L, U, T, X, or V.
    # We use a regular expression to match fare basis codes ("bc") that start with one of these letters.
    # All letters are lowercased as requested.
    fare_component = "/f bc=~^[lutxv]"

    components = [airline_component, stops_component, fare_component]

    # 2. Find the optimal ordering of the components.

    # Generate all possible orderings (permutations).
    permutations = list(itertools.permutations(components))

    # Join the components in each permutation with a semicolon to form valid code strings.
    code_strings = [";".join(p) for p in permutations]

    # Sort the list of code strings to find the one with the highest lexicographical value.
    # The key `str.lower` ensures the sort is case-insensitive, and `reverse=True` sorts from highest to lowest.
    code_strings.sort(key=str.lower, reverse=True)

    # The first item in the reverse-sorted list is our answer.
    final_code = code_strings[0]
    final_code_parts = final_code.split(';')

    # 3. Print the final result and a detailed breakdown.

    print("--- Final Extension Code ---")
    print(final_code)
    print("\n--- Breakdown of the Code ---")
    print(f"The code string '{final_code}' is constructed from three commands, ordered to maximize lexicographical value:")
    print(f"1. '{final_code_parts[0]}': This command ensures the flight is direct. The number '0' specifies zero layovers.")
    print(f"2. '{final_code_parts[1]}': This command specifies that the fare must be on Delta Air Lines ('dl'). Using 'f dl' is the shortest way to write this.")
    print(f"3. '{final_code_parts[2]}': This command filters for the required fare classes (L, U, T, X, V) using a regular expression, making it eligible for the companion certificate.")

# Execute the function to find and print the solution.
solve_ita_matrix_code()