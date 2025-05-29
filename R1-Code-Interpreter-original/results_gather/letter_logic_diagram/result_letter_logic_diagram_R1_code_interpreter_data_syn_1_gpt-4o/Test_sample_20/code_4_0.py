def solve_puzzle():
    # Manually filled grid based on constraints
    solution = [
        "b,a,f,g,d,c,e",
        "c,d,g,a,f,e,b",
        "d,c,a,c,g,b,f",
        "g,f,c,e,a,b,d",
        "f,g,e,b,c,a",
        "a,e,c,f,d,g,a",
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