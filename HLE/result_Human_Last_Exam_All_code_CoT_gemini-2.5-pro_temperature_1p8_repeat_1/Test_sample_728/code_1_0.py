import requests
import urllib.parse

def find_passing_films():
    """
    Finds which films from a predefined list pass the Bechdel test.
    """
    films = [
        ("a", "Igla", "1988"),
        ("b", "Ghost Town", "2008"),
        ("c", "Girls Will Be Girls", "2003"),
        ("d", "War Dogs", "2016"),
        ("e", "Slither", "2006"),
        ("f", "John Dies at the End", "2012"),
        ("g", "Man Who Knew Too Much, The", "1934"),
        ("h", "Ladies In Retirement", "1941"),
        ("i", "Internship, The", "2013"),
        ("j", "Tinker Bell and the Lost Treasure", "2009"),
    ]

    passing_films_letters = []
    base_url = "http://bechdeltest.com/api/v1/getMoviesByTitle"

    for letter, title, year in films:
        # Clean up title for better API search (e.g., "Movie, The" -> "The Movie")
        search_title = title
        if title.endswith(', The'):
            search_title = 'The ' + title[:-5]
            
        params = {'title': search_title}
        
        try:
            response = requests.get(base_url, params=params, timeout=10)
            response.raise_for_status()  # Raise an exception for bad status codes
            movies = response.json()
            
            for movie in movies:
                # The API returns year as a string, so we compare strings
                if str(movie.get('year')) == year:
                    # Rating 3 means it passes all 3 criteria
                    if movie.get('rating') == 3:
                        passing_films_letters.append(letter)
                    break # Found the correct movie for the year, move to the next film
                        
        except requests.exceptions.RequestException as e:
            # Silently ignore errors for this task to avoid clutter
            # print(f"Could not retrieve data for {title}: {e}")
            pass

    # Print the final result in the specified format
    result_string = ",".join(sorted(passing_films_letters))
    print(result_string)
    print(f"<<<{result_string}>>>")

if __name__ == "__main__":
    find_passing_films()