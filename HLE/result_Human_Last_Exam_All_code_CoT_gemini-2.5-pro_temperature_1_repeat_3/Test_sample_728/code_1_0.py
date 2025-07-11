import requests
import urllib.parse

def check_bechdel_test():
    """
    Checks a predefined list of films against the bechdeltest.com API
    to see if they pass all three conditions of the test.
    """
    # The list of films to evaluate
    films = {
        'a': ('Igla', '1988'),
        'b': ('Ghost Town', '2008'),
        'c': ('Girls Will Be Girls', '2003'),
        'd': ('War Dogs', '2016'),
        'e': ('Slither', '2006'),
        'f': ('John Dies at the End', '2012'),
        'g': ('The Man Who Knew Too Much', '1934'),
        'h': ('Ladies In Retirement', '1941'),
        'i': ('The Internship', '2013'),
        'j': ('Tinker Bell and the Lost Treasure', '2009')
    }

    passing_films_letters = []
    base_url = "http://bechdeltest.com/api/v1/getMoviesByTitle"

    print("Evaluating films based on the Bechdel Test criteria...")
    print("-" * 50)

    for letter, (title, year) in films.items():
        try:
            # Prepare the URL for the API request
            params = {'title': title}
            response = requests.get(base_url, params=params)
            response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)
            
            movies_data = response.json()
            
            found_movie = None
            # The API returns a list, so we find the correct movie by matching the year
            for movie in movies_data:
                if str(movie.get('year')) == year:
                    found_movie = movie
                    break
            
            # If no exact match, take the first result as a fallback
            if not found_movie and movies_data:
                found_movie = movies_data[0]
                
            if found_movie:
                rating = found_movie.get('rating', 0)
                if rating == 3:
                    status = "Pass"
                    passing_films_letters.append(letter)
                else:
                    status = f"Fail (Rating: {rating})"
                print(f"{letter}) {title} ({year}): {status}")
            else:
                print(f"{letter}) {title} ({year}): Not found in the database.")

        except requests.exceptions.RequestException as e:
            print(f"Could not retrieve data for '{title}': {e}")
        except ValueError:
            print(f"Could not parse response for '{title}'.")

    print("-" * 50)
    
    # Sort the letters alphabetically for a consistent result
    passing_films_letters.sort()
    
    # Format the final answer as a comma-separated string
    final_answer = ",".join(passing_films_letters)
    
    print("\nList of letters for films that pass all three conditions:")
    print(final_answer)
    
    # The final answer format required by the user
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    check_bechdel_test()