# 1. Define the components of the ITA Matrix extension code based on the requirements.

# The airline (Delta) and stop requirement (non-stop).
# 'f' specifies the marketing carrier, and can be lowercased. '/ nonstop' is also lowercased.
carrier_and_stops_code = "f DL / nonstop"

# The fare class requirement for the Platinum Delta Companion Certificate.
# 'f' specifies the booking codes. The fare classes L, U, T, X, and V must be uppercase.
eligible_fare_classes = ["L", "U", "T", "X", "V"]
fare_class_code = f"f {'|'.join(eligible_fare_classes)}"

# The separator used in ITA Matrix extension codes.
separator = " ; "

# 2. Combine the components to form the final, optimized extension code.
# This format was chosen because it is one of the shortest valid options and is
# lexicographically highest among them.
final_extension_code = f"{carrier_and_stops_code}{separator}{fare_class_code}"

# 3. Print the final result.
# The output displays the complete "equation" with all its parts (including the fare class "numbers").
print(final_extension_code)