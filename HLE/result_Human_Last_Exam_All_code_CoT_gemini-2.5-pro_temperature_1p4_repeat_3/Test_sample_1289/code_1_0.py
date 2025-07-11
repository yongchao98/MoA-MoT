# Step 1: Define the 7 intervals of the ascending Maqam Zanjaran.
# Jins Hijaz Kar on C -> [0.5, 1.5, 0.5]
# Linking interval F to G -> [1.0]
# Jins Rast on G -> [1.0, 0.75, 0.75]
ascending_intervals = [0.5, 1.5, 0.5, 1.0, 1.0, 0.75, 0.75]

# Step 2: Define the 4 intervals of the descending phrase.
# The descending notes are C', Bb, Ab, G, F.
# Interval C' to Bb -> 1.0
# Interval Bb to Ab -> 1.0
# Interval Ab to G -> 0.5
# Interval G to F -> 1.0
descending_intervals = [1.0, 1.0, 0.5, 1.0]

# Step 3: Combine the ascending and descending intervals into a single list.
all_intervals = ascending_intervals + descending_intervals

# Step 4: Format the list into the specified string format '{n1,n2,n3,...}'.
# First, convert all numbers in the list to their string representation.
interval_strings = [str(interval) for interval in all_intervals]

# Next, join these strings with a comma.
formatted_string = "{" + ",".join(interval_strings) + "}"

# Step 5: Print the final result.
print(formatted_string)

# The final answer in the required format.
print(f"<<<{formatted_string}>>>")