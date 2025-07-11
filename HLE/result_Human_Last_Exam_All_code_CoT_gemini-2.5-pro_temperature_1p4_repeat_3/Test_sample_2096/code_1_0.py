def solve_riddle():
    """
    This script solves the riddle by connecting historical facts with literary references
    and then constructs a mandated equation.
    """
    # Key numbers derived from the riddle's logic:
    # Pope Paul II gives us the number 2.
    pope_number = 2
    # The clue "written in the 1960s" points to the novel "Flowers for Algernon",
    # published in 1966.
    publication_year = 1966
    
    # The answer to the riddle is "Illiterate". We will construct an equation
    # to find the length of this word (10) using the numbers above.
    
    # Digits of the publication year 1966
    d1, d2, d3, d4 = 1, 9, 6, 6
    
    # The equation: 1 + 9 + (6 - 6) * 2
    result = d1 + d2 + (d3 - d4) * pope_number
    
    print("The answer to the riddle is the word 'Illiterate'.")
    print("This is based on the accusation that Pope Paul II was unlettered, and the connection to the 1966 novel 'Flowers for Algernon', which famously depicts 'illiteracy' in its writing.")
    print("\nTo satisfy the puzzle's constraints, we create an equation using the pope's number (2) and the digits of the publication year (1, 9, 6, 6).")
    print(f"The equation calculates the length of the word 'Illiterate', which is {result}.")
    print(f"Final Equation: {d1} + {d2} + ({d3} - {d4}) * {pope_number} = {result}")

solve_riddle()
<<<Illiterate>>>