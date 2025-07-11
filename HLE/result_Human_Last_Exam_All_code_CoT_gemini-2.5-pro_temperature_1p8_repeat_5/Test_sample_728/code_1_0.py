import requests
import urllib.parse

def check_bechdel_test():
    """
    Checks a predefined list of films against the Bechdel Test API
    and prints the films that pass.
    """
    films = [
        ('a', 'Igla', 1988),
        ('b', 'Ghost Town', 2008),
        ('c', 'Girls Will Be Girls', 2003),
        ('d', 'War Dogs', 2016),
        ('e', 'Slither', 2006),
        ('f', 'John Dies at the End', 2012),
        ('g', 'Man Who Knew Too Much, The', 1934),
        ('h', 'Ladies In Retirement', 1941),
        ('i', 'Internship, The', 2013),
        ('j', 'Tinker Bell and the Lost Treasure', 2009),
    ]

    base_url = "http://bechdeltest.com/api/v1/getMoviesByTitle"
    passing_films_letters = []

    print("Evaluating films based on the Bechdel Test...\n")

    for letter, title, year in films:
        # The API search is better if ", The" at the end is moved to the front.
        search_title = title
        if title.endswith(", The"):
            search_title = "The " + title[:-5]
        
        params = {'title': search_title}
        
        try:
            response = requests.get(base_url, params=params, timeout=10)
            response.raise_for_status()  # Will raise an HTTPError for bad responses (4xx or 5xx)
            movies = response.json()
            
            found_movie = None
            for movie in movies:
                if movie.get('year') == year:
                    found_movie = movie
                    break
            
            if found_movie:
                rating = found_movie.get('rating', -1)
                passes = "Passes" if rating == 3 else "Fails"
                print(f"{letter}) {title} ({year}) -> Rating: {rating} ({passes})")
                if rating == 3:
                    passing_films_letters.append(letter)
            else:
                print(f"{letter}) {title} ({year}) -> Could not find an exact match in the database.")

        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for '{title}': {e}")

    # Sort the letters alphabetically for a consistent output
    passing_films_letters.sort()
    
    # Final answer in the specified format
    final_answer = ",".join(passing_films_letters)
    print("\nFinal list of passing films (by letter):")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    check_bechdel_test()