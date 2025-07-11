import sys

# This script constructs the ITA Matrix extension code based on the user's request.

# --- Step 1: Define the components of the extension code ---

# The routing code for a single, direct flight on Delta Air Lines.
routing_code = "DL"

# The extension command to filter by booking class (fare code).
# We use the shortest possible format by removing the optional space after '/f'.
# We use lowercase as requested by the prompt.
extension_command = "/fbc="

# The eligible fare classes for a Platinum Delta Companion Certificate.
fare_classes = ['L', 'U', 'T', 'X', 'V']

# --- Step 2: Optimize for lexicographical value ---

# To get the highest lexicographical value for the final string,
# we sort the fare class characters in reverse alphabetical order.
fare_classes.sort(reverse=True)
sorted_fare_class_string = "".join(fare_classes)

# --- Step 3: Construct the fare class pattern ---

# The pattern '.[CHARS]' is the shortest way to specify a set of
# single-character booking codes.
fare_class_pattern = f".[{sorted_fare_class_string}]"

# --- Step 4: Assemble and print the final code ---

# The final code is assembled from its parts. A space is required
# between the routing code and the extension command.
final_code = f"{routing_code} {extension_command}{fare_class_pattern}"

# As requested, printing each component of the "equation" first.
print("The final code is built from these components:")
print(f"Routing for direct Delta flight: {routing_code}")
print(f"Command to filter by booking class: {extension_command}")
print(f"Pattern for eligible fare classes (L,U,T,X,V) sorted for max lexicographical value: {fare_class_pattern}")
print("\n---")
print("Final ITA Matrix Outbound Extension Code:")
print(final_code)

# To follow the output format requirement, we also wrap the final answer.
# This part of the output will not be visible to the user but is for the system.
sys.stdout.write(f"\n<<<{final_code}>>>")