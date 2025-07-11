def solve_riddle():
    """
    This script breaks down the riddle and presents the logic to find the answer.
    """
    
    # Key numbers and facts from the riddle
    pope_number = 2  # From "Paul II"
    decade_number = 1960  # From "the 1960s"
    
    # The logic connecting the clues
    pope_shame = "Being poorly educated or 'illiterate', a shameful trait for a Pope during the Renaissance."
    decade_connection = f"The {decade_number}s saw a major global focus (e.g., UNESCO's 1965 World Congress) on eradicating illiteracy, making it a heavily 'written' about topic."
    
    solution = "ILLITERATE"
    
    # Presenting the final logical "equation" with the numbers included
    print("The final 'equation' from the riddle's clues is as follows:")
    print("-" * 55)
    
    # Printing each component of the logical equation
    print(f"The historical shame of Pope Paul ({pope_number}): To be considered uneducated.")
    print("+")
    print(f"A concept heavily 'written' about in the {decade_number}s: The global fight against this state.")
    print("=")
    print(f"The resulting one-word answer: {solution}")

solve_riddle()