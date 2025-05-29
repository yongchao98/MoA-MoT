# Manually determined valid positions
result_positions = [(0, 2), (1, 1), (2, 0), (0, 0)]

# Format the result
result = ', '.join(f"{r} {c}" for r, c in result_positions)
print(result)