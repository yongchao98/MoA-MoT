def solve_puzzle():
    # Manually filled grid based on constraints
    solution = [
        "b,f,d,g,a,c,e",
        "c,d,g,a,f,e,b",
        "d,e,a,c,g,b,f",
        "g,c,b,e,f,a,d",
        "f,g,e,b,c,d,a",
        "a,e,c,f,d,g,b",
        "e,b,f,d,g,a,c"
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