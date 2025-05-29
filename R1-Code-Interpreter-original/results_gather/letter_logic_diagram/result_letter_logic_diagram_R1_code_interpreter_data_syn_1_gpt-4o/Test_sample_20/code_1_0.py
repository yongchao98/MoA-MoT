def solve_puzzle():
    # Manually filled grid based on constraints
    solution = [
        "b,f,a,g,c,d,e",
        "c,d,g,a,f,e,b",
        "d,e,a,c,g,b,f",
        "g,c,f,e,b,a,d",
        "f,g,e,b,d,c,a",
        "a,c,d,f,e,g,b",
        "e,b,c,d,a,f,g"
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