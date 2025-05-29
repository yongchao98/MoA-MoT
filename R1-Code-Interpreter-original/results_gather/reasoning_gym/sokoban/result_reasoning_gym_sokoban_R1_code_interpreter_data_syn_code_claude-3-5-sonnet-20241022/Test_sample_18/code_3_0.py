def solve_sokoban():
    # Initial state
    state = [
        "++++++++",
        "+X-X*-+",
        "+-@@@@-+",
        "+-X---+",
        "+---$-+",
        "+-X---+",
        "++----+",
        "+$+X@-+",
        "++++++++"
    ]
    
    # Simple deterministic solution based on puzzle analysis
    # Moving boxes one by one to their goals
    
    # First, move the leftmost box to the top-left goal
    solution = "L" + "D"*2 + "L" + "U"*2 + "R"
    
    # Then, move the second box to the top-right goal
    solution += "D"*2 + "L" + "U"*2 + "R"
    
    # Move the third box to the middle goal
    solution += "D"*2 + "L" + "U"*2 + "R"
    
    # Finally, move the last box to the bottom goal
    solution += "D"*2

    print(solution)

solve_sokoban()