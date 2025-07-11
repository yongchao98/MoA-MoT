import requests
import json

def format_title(title):
    """
    Reformats titles ending in ', The' to 'The ...' for better API matching.
    e.g., 'Internship, The' becomes 'The Internship'.
    """
    if title.endswith(', The'):
        return 'The ' + title[:-5]
    return title

def check_bechdel_test():
    """
    Checks a predefined list of films against the Bechdel Test API,
    prints the results for each, and then a final summary.
    """
    films = {
        'a': ('Igla', 1988),
        'b': ('Ghost Town', 2008),
        'c': ('Girls Will Be Girls', 2003),
        'd': ('War Dogs', 2016),
        'e': ('Slither', 2006),
        'f': ('John Dies at the End', 2012),
        'g': ('Man Who Knew Too Much, The', 1934),
        'h': ('Ladies In Retirement', 1941),
        'i': ('Internship, The', 2013),
        'j': ('Tinker Bell and the Lost Treasure', 2009)
    }

    base_url = "http://bechdeltest.com/api/v1/getMoviesByTitle"
    passing_films = []

    print("Checking films against the Bechdel Test database...")
    
    # Use a session for potential connection pooling
    with requests.Session() as session:
        for letter, (title, year) in films.items():
            formatted_title = format_title(title)
            print(f"Checking '{formatted_title}' ({year})...", end=' ')
            
            try:
                params = {'title': formatted_title}
                response = session.get(base_url, params=params, timeout=10)
                response.raise_for_status()  # Raise an exception for bad status codes
                
                results = response.json()
                
                film_found = False
                for movie in results:
                    # Find the exact movie by matching the year
                    if movie.get('year') == year:
                        rating = movie.get('rating')
                        if rating == 3:
                            print(f"Passed! (Rating: {rating})")
                            passing_films.append(letter)
                        else:
                            print(f"Failed. (Rating: {rating})")
                        film_found = True
                        break # Stop searching once the correct year is found
                
                if not film_found:
                    print(f"Failed. (Movie not found in database for year {year})")

            except requests.exceptions.RequestException as e:
                print(f"Error. (Could not connect to API: {e})")
            except json.JSONDecodeError:
                print("Error. (Could not parse API response)")

    passing_films.sort()
    result_string = ",".join(passing_films)
    
    print("\n---")
    print("The letters of the films that pass all three conditions of the Bechdel Test are:")
    print(f"<<<{result_string}>>>")


if __name__ == '__main__':
    check_bechdel_test()