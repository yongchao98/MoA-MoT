def solve_movie_query():
    """
    Identifies the film and provides details along with a
    numerical equation based on the film's award history.
    """
    film_title = "The Fugitive"
    character = "Dr. Charles Nichols"
    context = "He says 'Thank you' while getting on a Chicago bus, which ultimately leads to his capture by Dr. Richard Kimble."

    # To satisfy the prompt's requirement for an equation,
    # we'll use the film's Oscar statistics: 7 nominations and 1 win.
    oscar_nominations = 7
    oscar_wins = 1
    other_nominations = oscar_nominations - oscar_wins

    print(f"The film you are asking about is '{film_title}'.")
    print(f"The character is {character}, who seals his fate with this action.")
    print(context)
    print("\nAs requested, here is an equation related to the film's Oscar performance:")
    # The final print statement outputs each number in the equation.
    print(f"Nominations ({oscar_nominations}) - Unsuccessful Nominations ({other_nominations}) = Wins ({oscar_wins})")

solve_movie_query()