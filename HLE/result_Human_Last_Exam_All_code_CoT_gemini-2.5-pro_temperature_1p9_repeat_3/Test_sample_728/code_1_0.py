import requests
import sys

def check_bechdel_test():
    """
    Checks a predefined list of films against the Bechdel Test API
    and prints the films that pass.
    """
    # Note: If you are behind a proxy, you might need to configure it for the 'requests' library.
    # The 'requests' library needs to be installed (pip install requests).

    films = [
        {'letter': 'a', 'title': 'Igla', 'year': 1988},
        {'letter': 'b', 'title': 'Ghost Town', 'year': 2008},
        {'letter': 'c', 'title': 'Girls Will Be Girls', 'year': 2003},
        {'letter': 'd', 'title': 'War Dogs', 'year': 2016},
        {'letter': 'e', 'title': 'Slither', 'year': 2006},
        {'letter': 'f', 'title': 'John Dies at the End', 'year': 2012},
        {'letter': 'g', 'title': 'Man Who Knew Too Much, The', 'year': 1934},
        {'letter': 'h', 'title': 'Ladies In Retirement', 'year': 1941},
        {'letter': 'i', 'title': 'Internship, The', 'year': 2013},
        {'letter': 'j', 'title': 'Tinker Bell and the Lost Treasure', 'year': 2009},
    ]

    passing_films_letters = []
    base_url = 'http://bechdeltest.com/api/v1/getMoviesByTitle'

    def normalize_title(title):
        """Handle titles in 'Title, The' format."""
        if title.lower().endswith(", the"):
            return "The " + title[:-5]
        return title

    print("--- Checking Films Against Bechdeltest.com API ---")

    for film in films:
        letter, title, year = film['letter'], film['title'], film['year']
        print(f"\nChecking film: {letter}) {title} ({year})")

        # Use normalized title for better search results
        search_title = normalize_title(title)
        params = {'title': search_title}

        try:
            response = requests.get(base_url, params=params, timeout=10)
            response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
            
            movies_found = response.json()
            film_match = None
            
            # Find the correct film by matching title and year
            for movie_data in movies_found:
                # API title search can be fuzzy, so we check if our search title is IN the result title
                # Also handle cases like "Man Who Knew Too Much, The" (1934) vs. (1956)
                if (search_title.lower() in movie_data.get('title', '').lower() and 
                    movie_data.get('year') == year):
                    film_match = movie_data
                    break
            
            if film_match:
                rating = film_match.get('rating')
                if rating == 3:
                    print(f"> PASSED: Found with Bechdel rating of 3.")
                    passing_films_letters.append(letter)
                else:
                    print(f"> FAILED: Found with Bechdel rating of {rating}, but a rating of 3 is required.")
            else:
                print(f"> FAILED: Movie not found in the database with the specified year.")

        except requests.exceptions.RequestException as e:
            print(f"> ERROR: Could not connect to the API. {e}", file=sys.stderr)
        except ValueError: # handles JSON decoding errors
            print(f"> ERROR: Could not parse API response.", file=sys.stderr)

    print("\n--- Results ---")
    if passing_films_letters:
        # Final equation/result output as requested.
        print("The following films pass all three conditions of the Bechdel Test:")
        passing_string = ", ".join([f"Film '{f['title']}' ({f['letter']})" for f in films if f['letter'] in passing_films_letters])
        print(passing_string)

        final_answer = ",".join(sorted(passing_films_letters))
        print(f"\nFinal Answer (as requested format):")
        print(final_answer)
    else:
        print("No films on the list passed the Bechdel Test.")
    print("-------------")


if __name__ == '__main__':
    check_bechdel_test()