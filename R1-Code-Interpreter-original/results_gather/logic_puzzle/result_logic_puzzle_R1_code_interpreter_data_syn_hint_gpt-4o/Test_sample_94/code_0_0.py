selected_positions = [
    (0, 0), (0, 4), (0, 5),
    (1, 0), (1, 1),
    (2, 1), (2, 3),
    (3, 0), (3, 4),
    (4, 0), (4, 4),
    (5, 1), (5, 4)
]

# Format the output as required
output = ', '.join(f"{r} {c}" for r, c in selected_positions)
print(f"<<<{output}>>>")