def solve_movie_riddle():
    """
    This function identifies the film based on the user's query and
    prints the answer along with an equation as requested.
    """
    # Step 1: Define the information about the film.
    movie_title = "L.A. Confidential"
    release_year = 1997

    # Step 2: Print the answer to the user's trivia question.
    print(f"The Oscar-nominated film in which a character seals his fate by saying 'Thank you' while boarding a bus is '{movie_title}'.")
    print("-" * 20)

    # Step 3: Fulfill the requirement to output an equation with numbers.
    # The equation will be based on the film's release year.
    num1 = 1900
    num2 = 97
    result = release_year  # 1900 + 97 = 1997

    # The final code must output each number in the final equation.
    print("As requested, here is an equation featuring the film's release year:")
    print(f"{num1} + {num2} = {result}")

# Execute the function to provide the answer.
solve_movie_riddle()