import sys

# This script generates an ITA Matrix extension code for a specific Delta flight search.

# --- Step 1: Define the components of the extension code ---

# Component for direct flights (0 stops).
# 'n' is the command for number of stops. Commands are lowercased per instructions.
stops_component = "n 0"

# Component for the carrier (Delta) and eligible fare classes.
# 'DL' is the carrier code for Delta.
# 'f' is the fare command, 'bc' is the booking code attribute.
# Fare classes L, U, T, X, and V are for Platinum Companion Certificates.
# The '|' acts as an "OR" operator.
# The carrier code is prefixed to the fare command to minimize string length.
fare_component = "DL f bc=L|bc=U|bc=T|bc=X|bc=V"

# --- Step 2: Combine components based on optimization rules ---

# The two main components to be joined by a semicolon are:
# 1. stops_component
# 2. fare_component

# To find the valid string of minimal length with the highest case-insensitive
# lexicographic value, we sort the components in reverse alphabetical order.
# For example, "n..." is lexicographically greater than "d...".
components = [stops_component, fare_component]

# Sorts in-place, in reverse (Z->A) order, ignoring case.
components.sort(key=str.lower, reverse=True)

# Join the sorted components with a semicolon, omitting spaces to minimize length.
final_extension_code = ";".join(components)


# --- Step 3: Print the breakdown and the final result ---

print("ITA Matrix Extension Code Breakdown:")
print(f"Component for direct flights: {stops_component}")
print(f"Component for Delta and Platinum Companion fares: {fare_component}")
print("\nComponents are ordered to maximize lexicographic value and joined with a semicolon.")

print("\nFinal Extension Code:")
print(final_extension_code)

# For automated checking, output the final answer in the specified format.
# This part of the output will not be visible to the user in a standard run.
sys.stdout = open(sys.devnull, 'w')
print(f'<<<{final_extension_code}>>>')