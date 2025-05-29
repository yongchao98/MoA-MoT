# Initial weights of the penguins
weights = {
    "Louis": 11,
    "Bernard": 13,
    "Vincent": 11,
    "Gwen": 15,
    "James": 12  # Added James
}

# Remove Bernard
weights.pop("Bernard")

# Calculate the cumulated weight
cumulated_weight = sum(weights.values())

print(cumulated_weight)