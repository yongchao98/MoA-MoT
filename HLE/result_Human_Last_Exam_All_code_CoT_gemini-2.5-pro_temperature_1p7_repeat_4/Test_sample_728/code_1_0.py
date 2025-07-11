import requests
import urllib.parse

def find_passing_films():
    """
    Checks a predefined list of films against the Bechdel Test using the
    bechdeltest.com API and prints the letters of the films that pass.
    """
    # The list of films to check, identified by a letter.
    # Format: { 'letter': ('Title', Year) }
    films = {
        'a': ('Igla', 1988),
        'b': ('Ghost Town', 2008),
        'c': ('Girls Will Be Girls', 2003),
        'd': ('War Dogs', 2016),
        'e': ('Slither', 2006),
        'f': ('John Dies at the End', 2012),
        'g': ('The Man Who Knew Too Much', 1934), # Title adjusted for easier searching
        'h': ('Ladies In Retirement', 1941),
        'i': ('The Internship', 2013),       # Title adjusted for easier searching
        'j': ('Tinker Bell and the Lost Treasure', 2009)
    }

    passing_film_letters = []
    base_url = "http://bechdeltest.com/api/v1/getMoviesByTitle"
    print("Checking films against the Bechdel Test database...")

    for letter, (title, year) in films.items():
        params = {'title': title}
        try:
            # Make a request to the API
            response = requests.get(base_url, params=params, timeout=10)
            response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)
            movie_results = response.json()

            found_film = False
            # Look for the specific movie by matching title and year
            for movie in movie_results:
                if str(movie.get('year')) == str(year) and movie.get('title').lower() == title.lower():
                    if movie.get('rating') == 3:
                        passing_film_letters.append(letter)
                        print(f"-> PASSED: '{title}' ({year})")
                    else:
                        print(f"-> FAILED: '{title}' ({year}) - Rating: {movie.get('rating')}/3")
                    found_film = True
                    break
            
            if not found_film:
                 print(f"-> INFO: Could not find an exact match for '{title}' ({year}) in the database.")

        except requests.exceptions.RequestException as e:
            print(f"-> ERROR: Could not fetch data for '{title}'. Reason: {e}")

    # Sort the letters alphabetically for a consistent result
    passing_film_letters.sort()

    # Print the final result in the requested format
    result_string = ",".join(passing_film_letters)
    print("\n---")
    print("Final Answer: The letters corresponding to films that pass all three conditions are:")
    print(result_string)


if __name__ == '__main__':
    find_passing_films()