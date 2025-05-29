# Known positions of the "✅" marks
known_positions = [
    (1, 14), (2, 7), (3, 16), (4, 12), (5, 13), (6, 8), (7, 2), (7, 19),
    (8, 6), (9, 11), (10, 15), (11, 9), (12, 4), (13, 5), (14, 1), (14, 18),
    (15, 10), (16, 3), (19, 14), (20, 8)
]

# Check for a modular pattern
def check_modular_pattern(positions, mod_value):
    results = []
    for x, y in positions:
        if (y - x) % mod_value == 0:
            results.append((x, y))
    return results

# Try different mod values
for mod_value in range(1, 21):
    modular_positions = check_modular_pattern(known_positions, mod_value)
    if len(modular_positions) == len(known_positions):
        print(f"Modular Pattern Found with mod value {mod_value}: {modular_positions}")