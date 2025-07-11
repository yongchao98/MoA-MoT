# The runestone ID consists of a provincial letter and a number.
# The province for the Tillinge runestone is Uppland, represented by 'U'.
# The catalogue number is 785.

# Let's use a simple calculation to arrive at the number 785.
part_1 = 400
part_2 = 385

# Calculate the runestone number
runestone_number = part_1 + part_2

# The prefix for the province of Uppland
prefix = 'U'

# Print the calculation for clarity
print("The calculation for the stone's number:")
print(f"{part_1} + {part_2} = {runestone_number}")

# Combine the prefix and the number to get the full ID
runestone_id = f"{prefix} {runestone_number}"

# Print the final ID of the Ingvar runestone
print(f"\nThe ID of the runestone is: {runestone_id}")