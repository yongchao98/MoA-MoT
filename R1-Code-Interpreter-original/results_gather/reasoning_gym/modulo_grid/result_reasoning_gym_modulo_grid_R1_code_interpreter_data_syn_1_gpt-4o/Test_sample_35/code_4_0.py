# Define the positions of the "✅" symbols
positions = [
    (1, 14), (2, 8), (4, 11), (8, 3), (8, 17), 
    (9, 13), (12, 5), (12, 19), (13, 9), (14, 8), 
    (15, 2), (15, 16), (17, 14), (18, 8), (20, 12)
]

# Check for a pattern in the positions
def find_pattern(positions):
    # Check if there is a consistent pattern in the row and column indices
    for i, (row, col) in enumerate(positions):
        # Check if the row or column follows a specific sequence or rule
        # For example, check if they are part of an arithmetic or geometric sequence
        print(f"Index: {i}, Row: {row}, Col: {col}")

find_pattern(positions)