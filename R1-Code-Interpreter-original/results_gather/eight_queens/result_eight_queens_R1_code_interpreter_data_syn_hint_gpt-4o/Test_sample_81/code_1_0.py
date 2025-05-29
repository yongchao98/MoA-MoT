# Convert the list of tuples to the required format
queen_positions = [(0, 0), (1, 5), (2, 7), (3, 2), (4, 6), (5, 3), (6, 1), (7, 4)]
formatted_positions = ', '.join(f"{row} {col}" for row, col in queen_positions)
print(f"<<<{formatted_positions}>>>")