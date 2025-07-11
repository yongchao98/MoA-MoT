def find_movie_by_quote():
    """
    Identifies the movie based on the "Thank you" on a bus scene.
    """
    # Information about the film and the specific scene
    movie_title = "The Fugitive"
    year = 1993
    character_saying_line = "Dr. Charles Nichols"
    protagonist = "Dr. Richard Kimble"

    # The explanation of the scene
    explanation = (
        f"The Oscar-nominated film is '{movie_title}' from the year {year}.\n\n"
        f"In a key scene, the villain, {character_saying_line}, boards an 'L' train (often mistaken for a bus) in Chicago. "
        f"He says a simple 'Thank you' to the operator. The protagonist, {protagonist}, who is already on the train but hasn't "
        f"seen him yet, recognizes his voice. This recognition leads directly to the final confrontation."
    )

    # Print the final answer
    print(explanation)

# Execute the function to provide the answer
find_movie_by_quote()