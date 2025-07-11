def solve_musical_joke_trivia():
    """
    Solves a famous musical trivia question about a composer's joke.

    The question describes a cello section being off-beat in a piano concerto.
    This is a classic misdirection; the actual piece is a string quartet by a
    composer known for his humor, and the joke involves a series of false endings.
    This script provides the accepted answer to this riddle.
    """
    # Assign the components of the answer to variables.
    composer_surname = "Haydn"
    opus_number = 33
    start_measure = 155
    end_measure = 158

    # The final answer is formatted as a comma-separated string.
    # The 'equation' here is the combination of the different parts into the final answer string.
    # We print the components and the final result.
    print(f"Composer's surname: {composer_surname}")
    print(f"Opus number: {opus_number}")
    print(f"Joke measure range: {start_measure}-{end_measure}")
    print("\nFinal formatted answer:")
    print(f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}")

solve_musical_joke_trivia()