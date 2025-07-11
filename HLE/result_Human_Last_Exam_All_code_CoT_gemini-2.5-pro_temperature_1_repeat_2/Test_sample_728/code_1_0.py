import requests
import urllib.parse

def find_passing_films():
    """
    Evaluates a predefined list of films against the Bechdel Test using an online API.
    Prints a comma-separated list of the letters corresponding to films that pass.
    """
    # A dictionary mapping letters to film titles and years
    films = {
        "a": ("Igla", 1988),
        "b": ("Ghost Town", 2008),
        "c": ("Girls Will Be Girls", 2003),
        "d": ("War Dogs", 2016),
        "e": ("Slither", 2006),
        "f": ("John Dies at the End", 2012),
        "g": ("Man Who Knew Too Much, The", 1934),
        "h": ("Ladies In Retirement", 1941),
        "i": ("Internship, The", 2013),
        "j": ("Tinker Bell and the Lost Treasure", 2009),
    }

    passing_films = []
    base_url = "http://bechdeltest.com/api/v1/getMoviesByTitle?title="

    # Iterate through each film to check its status
    for letter, (title, year) in films.items():
        # The API search can be particular about articles like "The"
        search_title = title
        if title.endswith(", The"):
            search_title = "The " + title.replace(", The", "")
        
        # Properly encode the title for the URL
        encoded_title = urllib.parse.quote_plus(search_title)
        api_url = base_url + encoded_title
        
        try:
            # Make the API request with a timeout
            response = requests.get(api_url, timeout=10)
            response.raise_for_status()  # Raise an error for bad responses (4xx or 5xx)
            movie_results = response.json()
            
            # Find the specific movie by matching its release year
            for movie in movie_results:
                if 'year' in movie and movie['year'] == str(year):
                    # A rating of 3 signifies that all Bechdel Test criteria are met
                    if 'rating' in movie and movie['rating'] == 3:
                        passing_films.append(letter)
                    break # Move to the next film once the correct one is found

        except requests.exceptions.RequestException as e:
            # Handle cases where the API call fails
            print(f"Warning: Could not fetch data for '{title}'. Error: {e}")
        except ValueError:
            # Handle cases where the JSON response is invalid
            print(f"Warning: Could not decode JSON for '{title}'.")

    # Sort the letters alphabetically and format the final string
    passing_films.sort()
    final_answer = ",".join(passing_films)
    
    print(final_answer)

if __name__ == "__main__":
    find_passing_films()