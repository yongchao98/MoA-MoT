def solve_trivia():
    """
    This function answers a historical question about Rome's Piazza della Rotonda
    and demonstrates a simple equation as requested.
    """
    # The year the feature was removed.
    year_removed = 1950
    # An arbitrary number for our equation. Let's use the number of rails for a tram track.
    number_of_rails = 2
    
    # An arbitrary calculation to fulfill the request for an equation.
    result = year_removed - number_of_rails

    # The answer to the user's question.
    feature_description = "A circular tram line that went around the square."
    
    print("The unique feature in Piazza della Rotonda, removed around 1950, was:")
    print(feature_description)
    print("\nTo fulfill the request for an equation, here is one using the year of removal:")
    print(f"Equation: {year_removed} - {number_of_rails} = {result}")

solve_trivia()