# Define the items on the table
items = {
    "black": {"mug": 1, "puzzle": 3, "dog_leash": 2, "keychain": 1},
    "pink": {"keychain": 1, "puzzle": 3},
    "mauve": {"keychain": 3, "mug": 3, "puzzle": 1}
}

# Count the number of mauve items
mauve_items = sum(items["mauve"].values())

# Print the result
print(mauve_items)