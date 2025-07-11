import requests

def evaluate_bechdel_test():
    """
    Evaluates a predefined list of films against the Bechdel Test
    by querying the bechdeltest.com API.
    """
    # 1. Define the list of films with their titles and years.
    # Titles have been slightly normalized for better API matching.
    films = {
        'a': ('Igla', 1988),
        'b': ('Ghost Town', 2008),
        'c': ('Girls Will Be Girls', 2003),
        'd': ('War Dogs', 2016),
        'e': ('Slither', 2006),
        'f': ('John Dies at the End', 2012),
        'g': ('The Man Who Knew Too Much', 1934),
        'h': ('Ladies In Retirement', 1941),
        'i': ('The Internship', 2013),
        'j': ('Tinker Bell and the Lost Treasure', 2009)
    }

    passed_letters = []
    api_url = "https://bechdeltest.com/api/v1/getMoviesByTitle"

    print("Checking Bechdel Test ratings for each film...")
    print("-" * 50)

    # 3. For each film, query the API.
    for letter, (title, year) in sorted(films.items()):
        try:
            response = requests.get(api_url, params={'title': title}, timeout=15)
            response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
            
            search_results = response.json()
            film_found_and_checked = False
            
            # 4. Find the specific film by matching the year and evaluate its rating.
            for movie_data in search_results:
                if str(movie_data.get('year')) == str(year):
                    rating = movie_data.get('rating', 0)
                    if rating == 3:
                        passed_letters.append(letter)
                        print(f"PASS: '{title}' ({year}) passes all 3 tests.")
                    else:
                        print(f"FAIL: '{title}' ({year}) does not pass all 3 tests (Rating: {rating}).")
                    film_found_and_checked = True
                    break  # Found the correct film, so we can stop searching for this title
            
            if not film_found_and_checked:
                print(f"FAIL: '{title}' ({year}) could not be found with a matching year in the database.")

        except requests.exceptions.RequestException as e:
            print(f"ERROR: Could not retrieve data for '{title}'. Reason: {e}")
        except ValueError:
            print(f"ERROR: Could not parse the API response for '{title}'.")

    # 5. Compile and print the final result.
    final_answer = ",".join(passed_letters)
    
    print("-" * 50)
    print("The letters for the films which pass all three conditions are:")
    print(final_answer)


if __name__ == '__main__':
    evaluate_bechdel_test()
<<<c,e,f,g,h,j>>>