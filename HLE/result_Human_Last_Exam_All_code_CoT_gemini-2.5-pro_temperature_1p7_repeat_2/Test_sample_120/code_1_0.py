def generate_ita_matrix_code():
    """
    This script determines the ITA Matrix extension code for finding direct Delta
    flights eligible for Platinum Companion Certificates, applying specific
    formatting and optimization rules.
    """

    # 1. Define the individual extension code components.
    # Code for direct flights (0 stops). 's' is the alias for 'maxstops'.
    direct_flight_code = "s 0"
    # Code for Platinum Companion Certificate fare classes (L, U, T, X, V).
    # 'f' is the fare code, 'bc' is the booking class. Lowercase is used as requested.
    fare_class_code = "f bc=l|u|t|x|v"

    components = [direct_flight_code, fare_class_code]

    # 2. Optimize the code string.
    # To find the string of minimal length with the highest case-insensitive
    # lexicographic value, we sort the components in reverse lexicographical order.
    # Since 's' > 'f', the component 's 0' will come first.
    sorted_components = sorted(components, key=str.lower, reverse=True)

    # 3. Join the components with a semicolon to form the final code.
    final_code = ";".join(sorted_components)

    # 4. Print the components and the final result as requested.
    print("ITA Matrix Extension Code Components:")
    print(f"Direct Flight Rule: {sorted_components[0]}")
    print(f"Fare Class Rule: {sorted_components[1]}")
    print("\nFinal Optimized Extension Code:")
    print(final_code)

    return final_code

# Execute the function and capture the final answer for the required format.
final_answer = generate_ita_matrix_code()

# The final answer in the required format
# <<<s 0;f bc=l|u|t|x|v>>>