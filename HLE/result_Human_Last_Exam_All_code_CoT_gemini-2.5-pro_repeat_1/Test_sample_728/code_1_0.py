import json

def solve_bechdel_test():
    """
    Evaluates a list of films against the Bechdel Test and identifies which ones pass.
    A film passes if it has a rating of 3.
    """
    # Data representing the films and their Bechdel test scores.
    # A rating of 3 means it passes all three criteria.
    # This data is based on information from bechdeltest.com.
    films_data = [
        {"letter": "a", "title": "Igla", "year": 1988, "rating": 1},
        {"letter": "b", "title": "Ghost Town", "year": 2008, "rating": 1},
        {"letter": "c", "title": "Girls Will Be Girls", "year": 2003, "rating": 3},
        {"letter": "d", "title": "War Dogs", "year": 2016, "rating": 1},
        {"letter": "e", "title": "Slither", "year": 2006, "rating": 3},
        {"letter": "f", "title": "John Dies at the End", "year": 2012, "rating": 1},
        {"letter": "g", "title": "Man Who Knew Too Much, The", "year": 1934, "rating": 1},
        {"letter": "h", "title": "Ladies In Retirement", "year": 1941, "rating": 3},
        {"letter": "i", "title": "Internship, The", "year": 2013, "rating": 2},
        {"letter": "j", "title": "Tinker Bell and the Lost Treasure", "year": 2009, "rating": 3}
    ]

    passing_films_letters = []
    
    print("Evaluating films based on the Bechdel Test:")
    print("-------------------------------------------")
    print("A film passes if its rating is 3.")
    print("1. Two named women characters.")
    print("2. They talk to each other.")
    print("3. About something besides a man.\n")

    for film in films_data:
        passes = film["rating"] == 3
        status = "Passes" if passes else "Fails"
        print(f"{film['letter']}) {film['title']} ({film['year']}) - Rating: {film['rating']} -> {status}")
        if passes:
            passing_films_letters.append(film["letter"])
            
    # Format the final answer as a comma-separated string
    final_answer = ",".join(passing_films_letters)
    
    print("\nFilms that pass all three conditions correspond to the letters:")
    print(final_answer)
    print("<<<{}>>>".format(final_answer))

solve_bechdel_test()