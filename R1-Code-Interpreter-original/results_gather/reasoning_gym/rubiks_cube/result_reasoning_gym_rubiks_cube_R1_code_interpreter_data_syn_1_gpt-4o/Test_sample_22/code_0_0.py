# Solution in Singmaster notation
solution = [
    # Solve the white cross
    "F R U R' U' F'",

    # Solve the white corners
    "U R U' L' U R' U' L",

    # Solve the middle layer edges
    "U R U' R' U' F' U F",

    # Solve the yellow cross
    "F R U R' U' F'",

    # Position the yellow edges
    "R U R' U R U2 R'",

    # Position the yellow corners
    "U R U' L' U R' U' L",

    # Orient the yellow corners
    "R U R' U R U2 R'"
]

# Print the solution
print(" ".join(solution))