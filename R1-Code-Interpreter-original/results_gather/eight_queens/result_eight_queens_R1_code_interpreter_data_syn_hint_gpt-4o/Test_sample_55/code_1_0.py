def format_queen_positions(positions):
    formatted_positions = ', '.join(f"{row} {col}" for row, col in positions)
    return formatted_positions

queen_positions = [(0, 0), (1, 4), (2, 7), (3, 5), (4, 2), (5, 6), (6, 1), (7, 3)]
formatted_positions = format_queen_positions(queen_positions)
print(formatted_positions)