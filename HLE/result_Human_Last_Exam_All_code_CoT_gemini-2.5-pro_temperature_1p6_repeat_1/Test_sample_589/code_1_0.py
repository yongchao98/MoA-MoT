import textwrap

def find_movie_scene():
    """
    Identifies and describes a specific movie scene based on a detailed query.
    """
    film_title = "Now and Then"
    release_year = 1995
    studio = "New Line Cinema"
    character = "Chrissy Dewitt"
    scene_description = (
        "In the 1995 film 'Now and Then', there is a scene where the four main "
        "characters confront the local bullies, the Wormer brothers. During the "
        "standoff, the character Chrissy Dewitt (who is 12 years old) panics "
        "while on her pink Schwinn bicycle and pedals directly toward the "
        "Wormer brothers' oncoming car. This forces the car to swerve off the "
        "road and into a ditch."
    )

    print(f"Movie Title: {film_title}")
    print(f"Release Year: {release_year}")
    print(f"Studio: {studio}")
    print("\nScene Details:")
    # Wrap the text for better readability in the console
    print(textwrap.fill(scene_description, width=80))

if __name__ == "__main__":
    find_movie_scene()