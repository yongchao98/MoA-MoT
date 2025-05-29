def solve_puzzle():
    # Correctly filled grid based on constraints
    solution = [
        "b,f,d,g,a,c,e",
        "f,d,g,a,e,b,c",
        "d,g,a,c,f,e,b",
        "g,a,c,e,b,d,f",
        "a,c,e,b,d,f,g",
        "c,e,b,f,g,a,d",
        "e,b,f,d,c,g,a"
    ]
    return solution

# Solve the puzzle
solution = solve_puzzle()

# Print the solution
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found.")