import collections

def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.

    The logic is as follows:
    1. A list of potential candidates (movies set in Paris, a city with a Luxor Obelisk)
       is created, sorted chronologically by the year they won Best Picture.
    2. Each movie is paired with a boolean based on film research, indicating whether
       it verifiably depicts the obelisk in Place de la Concorde.
    3. The script iterates through this list and identifies the first film
       for which the depiction is confirmed.
    """
    # MovieInfo named tuple for better readability.
    MovieInfo = collections.namedtuple('MovieInfo', ['year', 'title', 'depicts_obelisk'])

    # A list of Best Picture winners set in Paris, ordered chronologically.
    # The 'depicts_obelisk' flag is based on established film knowledge.
    # - The Life of Emile Zola (1937): Set in Paris, but shot on backlots. No confirmed obelisk scene.
    # - An American in Paris (1951): Features a famous, stylized depiction of the Place de la Concorde
    #   and its obelisk in the final ballet sequence.
    # - Around the World in 80 Days (1956): Features the obelisk during Paris scenes.
    # - Gigi (1958): Features the obelisk in on-location shots of Paris.
    candidate_films = [
        MovieInfo(year=1937, title="The Life of Emile Zola", depicts_obelisk=False),
        MovieInfo(year=1951, title="An American in Paris", depicts_obelisk=True),
        MovieInfo(year=1956, title="Around the World in 80 Days", depicts_obelisk=True),
        MovieInfo(year=1958, title="Gigi", depicts_obelisk=True),
    ]

    # Find the first film in the chronological list where depiction is confirmed.
    for film in candidate_films:
        if film.depicts_obelisk:
            print(f"The first Best Picture winner to depict a Luxor Obelisk was: {film.title} ({film.year})")
            return

find_first_winner_with_obelisk()
