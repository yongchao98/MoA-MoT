# Known positions of the "✅" marks
known_positions = [
    (1, 14), (2, 7), (3, 16), (4, 12), (5, 13), (6, 8), (7, 2), (7, 19),
    (8, 6), (9, 11), (10, 15), (11, 9), (12, 4), (13, 5), (14, 1), (14, 18),
    (15, 10), (16, 3), (19, 14), (20, 8)
]

# Analyze the differences between x and y coordinates
differences = [(x, y, y - x) for x, y in known_positions]

# Check for any consistent pattern in differences
def analyze_differences(differences):
    positive_diffs = [diff for _, _, diff in differences if diff > 0]
    negative_diffs = [diff for _, _, diff in differences if diff < 0]
    zero_diffs = [diff for _, _, diff in differences if diff == 0]
    return positive_diffs, negative_diffs, zero_diffs

positive_diffs, negative_diffs, zero_diffs = analyze_differences(differences)

# Print the analysis
print("Positive Differences:", positive_diffs)
print("Negative Differences:", negative_diffs)
print("Zero Differences:", zero_diffs)