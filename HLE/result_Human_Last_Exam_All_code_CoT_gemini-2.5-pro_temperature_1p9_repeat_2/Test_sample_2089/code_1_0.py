def solve_trivia():
    """
    This function provides the answer to the film trivia question and
    displays an equation as requested.
    """
    
    # Information about the film
    film_title = "The Godfather Part II"
    character_name = "Fredo Corleone"
    
    # Numbers for the equation
    release_year_part_2 = 1974
    release_year_part_1 = 1972
    resultant_part_number = 2
    
    # Print the answer to the user's question
    print(f"The Oscar-nominated film is '{film_title}'.")
    print(f"In the film, the character {character_name} seals his fate by saying \"Thank you\" while boarding a bus.")
    
    # Print the equation as requested
    print("\nHere is an equation using numbers related to the film series:")
    print(f"{release_year_part_2} - {release_year_part_1} = {resultant_part_number}")

# Execute the function
solve_trivia()