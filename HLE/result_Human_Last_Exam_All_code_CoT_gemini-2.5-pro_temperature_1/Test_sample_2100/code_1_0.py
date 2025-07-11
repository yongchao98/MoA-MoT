# Plan:
# 1. The clues provided point to Nikolai Rezanov, the protagonist of Andrei Voznesensky's rock opera "Juno and Avos".
# 2. The task is to provide his last name.
# 3. To meet the special requirement of an "equation with numbers", I will construct the name using the ASCII values of its letters.

# Define the ASCII values for each character in the last name "Rezanov".
char_codes = {
    'R': 82,
    'e': 101,
    'z': 122,
    'a': 97,
    'n': 110,
    'o': 111,
    'v': 118
}

# Construct the equation string for printing.
equation_parts = [f"chr({code})" for code in char_codes.values()]
equation_str = " + ".join(equation_parts)

print("The 'equation' to form the name from the character codes of its letters:")
print(equation_str)

# Build the final name by converting each code back to a character.
last_name = ""
for code in char_codes.values():
    last_name += chr(code)

print("\nThe final name is:")
print(last_name)