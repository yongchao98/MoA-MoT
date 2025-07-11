def evaluate_films():
    """
    This function evaluates a predefined list of films against the Bechdel Test
    and prints the letters of the films that pass.
    """

    # Film data including a 'passes_bechdel_test' flag.
    # This data is based on information from bechdeltest.com.
    film_data = {
        'a': {"title": "Igla (1988)", "passes_bechdel_test": False},
        'b': {"title": "Ghost Town (2008)", "passes_bechdel_test": False},
        'c': {"title": "Girls Will Be Girls (2003)", "passes_bechdel_test": True},
        'd': {"title": "War Dogs (2016)", "passes_bechdel_test": False},
        'e': {"title": "Slither (2006)", "passes_bechdel_test": True},
        'f': {"title": "John Dies at the End (2012)", "passes_bechdel_test": False},
        'g': {"title": "The Man Who Knew Too Much (1934)", "passes_bechdel_test": True},
        'h': {"title": "Ladies In Retirement (1941)", "passes_bechdel_test": True},
        'i': {"title": "The Internship (2013)", "passes_bechdel_test": False},
        'j': {"title": "Tinker Bell and the Lost Treasure (2009)", "passes_bechdel_test": True},
    }

    passing_letters = []
    print("Evaluating films for Bechdel Test compliance:")
    for letter, data in film_data.items():
        if data["passes_bechdel_test"]:
            passing_letters.append(letter)
            print(f"- {data['title']}: Passes")
        else:
            print(f"- {data['title']}: Fails")
    
    # Format the final result as a comma-separated list of letters
    result = ",".join(passing_letters)
    
    print("\nAnswer:")
    print(result)

evaluate_films()