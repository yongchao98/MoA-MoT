# Define the positions of the "✅" symbols
positions = [
    (1, 14), (2, 8), (4, 11), (8, 3), (8, 17), 
    (9, 13), (12, 5), (12, 19), (13, 9), (14, 8), 
    (15, 2), (15, 16), (17, 14), (18, 8), (20, 12)
]

# Check for a pattern in the positions
def find_pattern(positions):
    for row, col in positions:
        # Check if there is a consistent pattern
        # For example, check if the sum or difference of row and col is constant
        # Or if there is a consistent modulo operation
        # This is a placeholder for the actual pattern check
        print(f"Row: {row}, Col: {col}, Sum: {row + col}, Diff: {row - col}")

find_pattern(positions)