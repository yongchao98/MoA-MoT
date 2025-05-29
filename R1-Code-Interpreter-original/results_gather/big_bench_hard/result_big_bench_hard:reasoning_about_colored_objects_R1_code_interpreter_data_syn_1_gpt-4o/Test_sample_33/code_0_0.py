# Initial items on the table
items = {
    "green_bracelets": 3,
    "teal_dog_leash": 1,
    "green_dog_leash": 1,
    "green_paperclips": 3
}

# Remove teal items
items["teal_dog_leash"] = 0

# Count remaining paperclips
remaining_paperclips = items["green_paperclips"]

print(remaining_paperclips)