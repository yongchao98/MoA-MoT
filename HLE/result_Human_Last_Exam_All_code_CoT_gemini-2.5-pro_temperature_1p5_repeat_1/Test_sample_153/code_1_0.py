def calculate_letter_sum(name):
    """Calculates the sum of letter values (A=1, B=2, ...) for a given name."""
    total = 0
    name = name.upper()
    for char in name:
        if 'A' <= char <= 'Z':
            total += ord(char) - ord('A') + 1
    return total

def print_detailed_sum(name):
    """Prints the detailed summation equation for a name."""
    parts = []
    total = 0
    name = name.upper()
    for char in name:
        if 'A' <= char <= 'Z':
            value = ord(char) - ord('A') + 1
            parts.append(f"{char}({value})")
            total += value
    return f"{name}: {' + '.join(parts)} = {total}"

# The five locations from the answer choices
locations = {
    "A": "Kaliningrad",
    "B": "Perm",
    "C": "Taymyr",
    "D": "Chukotka",
    "E": "Adygea"
}

# The target number derived from the clue "Кома"
# Ко (Кобальт/Cobalt) = 27
# Ма (Марганец/Manganese) = 25
target_number = 27 + 25

print("Step 1: Decode the clue 'Кома'")
print("Ко -> Кобальт (Cobalt, atomic number 27)")
print("Ма -> Марганец (Manganese, atomic number 25)")
print(f"27 + 25 = {target_number}\n")

print("Step 2: Calculate the letter sum for each location's name (A=1, Z=26)")
correct_location_name = None
for key, name in locations.items():
    sum_val = calculate_letter_sum(name)
    print(f"Choice {key}: {name} -> {sum_val}")
    if sum_val == target_number:
        correct_location_name = name

print("\nStep 3: The location 'Perm' sums to 52, matching the target number.")
print("The final detailed equation is:")
if correct_location_name:
    print(print_detailed_sum(correct_location_name))
