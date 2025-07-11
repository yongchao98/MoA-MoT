import sys

def solve_ita_matrix_code():
    """
    This script determines the ITA Matrix extension code based on a set of rules.
    It constructs the code step-by-step and prints the logic and the final result.
    """

    # --- Step 1: Define the components based on certificate rules ---

    # Rule: Flight must be operated by Delta.
    # The standard code is O:DL. Per the prompt, we use lowercase where possible.
    operating_carrier_component = "o:dl"

    # Rule: Eligible fare classes are L, U, T, X, V.
    fare_classes = ["L", "U", "T", "X", "V"]


    # --- Step 2: Apply optimization and sorting rules ---

    # Rule: To maximize the lexicographical value, the fare classes must be
    # sorted in reverse alphabetical order.
    fare_classes.sort(reverse=True)
    
    # Construct the fare class component string.
    fare_class_component = "f " + "|".join(fare_classes)

    # The two components are 'o:dl' and 'f X|V|U|T|L'.
    # Rule: To maximize lexicographical value, the final string must start with
    # the component that is alphabetically higher ('o' > 'f').
    # Therefore, the operating carrier code comes first.
    first_part = operating_carrier_component
    second_part = fare_class_component
    separator = " ; "
    
    final_code = first_part + separator + second_part

    # --- Step 3: Print the explanation and the final answer ---

    print("To construct the ITA Matrix extension code, we follow these steps:")
    
    print("\n1. Define Constraints:")
    print("   - The flight must be operated by Delta, not just marketed by them.")
    print("   - The fare class must be one of L, U, T, X, or V.")
    
    print("\n2. Build Code Components:")
    print(f"   - Operating Carrier Code (lowercased): '{first_part}'")
    
    # We show the list of fare classes and how they are sorted for the final code.
    print(f"   - Eligible Fare Classes: {['L', 'U', 'T', 'X', 'V']}")
    print(f"   - Fare Classes sorted for max lexicographical value: {fare_classes}")
    print(f"   - Fare Class Code: '{second_part}'")
    
    print("\n3. Assemble Final Code:")
    print("   - The components are ordered to produce the highest lexicographical value ('o' comes after 'f').")
    print(f"   - Final Code Construction: '{first_part}' + '{separator}' + '{second_part}'")
    
    print("\n----------------------------------")
    print("Final ITA Matrix Extension Code:")
    print(final_code)
    print("----------------------------------")

    # Output the answer in the required format for the system.
    sys.stdout.write(f"\n<<<{final_code}>>>")

solve_ita_matrix_code()