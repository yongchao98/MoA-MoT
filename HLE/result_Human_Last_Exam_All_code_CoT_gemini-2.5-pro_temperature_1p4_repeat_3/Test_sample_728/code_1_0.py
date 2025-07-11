import collections

def solve_bechdel_test():
    """
    This function evaluates a predefined list of films against the Bechdel Test
    and prints the letters corresponding to the films that pass.
    """
    # Film data based on research. A film is marked 'True' if it passes all three Bechdel test criteria.
    # 1. Igla (1988): Fails. No two named women talk to each other.
    # 2. Ghost Town (2008): Passes. Gwen and Marjorie talk about the mummy.
    # 3. Girls Will Be Girls (2003): Passes. The main characters (Evie, Coco) talk about their careers and lives.
    # 4. War Dogs (2016): Fails. No two named women talk to each other.
    # 5. Slither (2006): Passes. Starla and Kylie talk about the alien parasite.
    # 6. John Dies at the End (2012): Fails. No two named women talk to each other.
    # 7. The Man Who Knew Too Much (1934): Fails. The two named women (Jill and Betty) only talk about Betty's father.
    # 8. Ladies In Retirement (1941): Passes. The women talk about the house, money, and their past.
    # 9. The Internship (2013): Fails. No two named women talk to each other.
    # 10. Tinker Bell and the Lost Treasure (2009): Passes. The female fairies discuss the moonstone and their quest.

    films = [
        {'letter': 'a', 'title': 'Igla', 'year': 1988, 'passes': False},
        {'letter': 'b', 'title': 'Ghost Town', 'year': 2008, 'passes': True},
        {'letter': 'c', 'title': 'Girls Will Be Girls', 'year': 2003, 'passes': True},
        {'letter': 'd', 'title': 'War Dogs', 'year': 2016, 'passes': False},
        {'letter': 'e', 'title': 'Slither', 'year': 2006, 'passes': True},
        {'letter': 'f', 'title': 'John Dies at the End', 'year': 2012, 'passes': False},
        {'letter': 'g', 'title': 'Man Who Knew Too Much, The', 'year': 1934, 'passes': False},
        {'letter': 'h', 'title': 'Ladies In Retirement', 'year': 1941, 'passes': True},
        {'letter': 'i', 'title': 'Internship, The', 'year': 2013, 'passes': False},
        {'letter': 'j', 'title': 'Tinker Bell and the Lost Treasure', 'year': 2009, 'passes': True}
    ]

    # Collect the letters of the films that pass the test
    passing_films_letters = []
    for film in films:
        if film['passes']:
            passing_films_letters.append(film['letter'])
    
    # Format the final answer as a comma-separated string
    answer = ",".join(passing_films_letters)
    
    print(answer)

solve_bechdel_test()