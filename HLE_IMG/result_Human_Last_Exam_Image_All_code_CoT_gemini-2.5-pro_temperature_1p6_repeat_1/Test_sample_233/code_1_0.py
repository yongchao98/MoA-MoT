# The runestone ID consists of a catalog prefix and a number.
# The catalog is Södermanlands runinskrifter, abbreviated as "Sö".
catalog_prefix = "Sö"

# The number of the runestone is 173.
# Let's construct this number from its digits to form an "equation".
hundreds = 1
tens = 7
ones = 3

# The equation to form the number
number = hundreds * 100 + tens * 10 + ones * 1

# The final ID is the concatenation of the prefix and the number.
runestone_id = f"{catalog_prefix} {number}"

# Print the "equation" as requested
print(f"The ID number is derived from the following calculation:")
print(f"{hundreds} * 100 + {tens} * 10 + {ones} * 1 = {number}")

# Print the final answer
print(f"\nThe ID of the Ingvar runestone is: {runestone_id}")
