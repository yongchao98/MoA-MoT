# Expanded solution for the given Rubik's Cube state
expanded_solution = [
    # Solve the White Cross
    "F R U R' U' F'",
    # Solve the White Corners
    "U R U' L' U R' U' L",
    # Solve the Middle Layer Edges
    "U R U' R' U' F' U F",
    # Solve the Yellow Cross
    "F R U R' U' F'",
    # Position the Yellow Edges
    "R U R' U R U U R'",
    # Position the Yellow Corners
    "U R U' L' U R' U' L",
    # Orient the Yellow Corners
    "R U R' U R U U R'"
]

# Convert the expanded solution into a single string with each move separated by a space
expanded_solution_string = ' '.join(expanded_solution)

# Print the expanded solution in Singmaster notation
print(expanded_solution_string)