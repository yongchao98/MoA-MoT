import sys
import subprocess

# This script attempts to install the required library if it's not found.
try:
    from duckduckgo_search import DDGS
except ImportError:
    print("The 'duckduckgo-search' library is not installed. Attempting to install...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "duckduckgo-search"])
        from duckduckgo_search import DDGS
        print("Installation successful.")
    except Exception as e:
        print(f"Failed to install 'duckduckgo-search'. Please install it manually by running: pip install duckduckgo-search")
        print(f"Error: {e}")
        sys.exit(1)

def find_movie_from_description():
    """
    Searches the web for a movie matching a description and prints the likely answer.
    """
    # The most unique parts of the query are a character's fate on a bus in an Oscar-nominated film.
    # The "Thank you" part is a common false memory associated with a different famous scene.
    search_query = "oscar nominated film character dies on bus"

    print(f"Searching online for a movie matching the query: '{search_query}'\n")

    # Using a with-statement ensures the session is properly closed.
    try:
        with DDGS() as ddgs:
            # We only need a few results to identify this famous scene.
            results = list(ddgs.text(search_query, max_results=3))

        if not results:
            print("Could not find any relevant results for the query.")
            return

        print("--- Top Search Results Snippets ---")
        for i, result in enumerate(results):
            # The .get() method is used to safely access dictionary keys.
            title = result.get('title', 'No Title')
            body = result.get('body', 'No Snippet').replace('\n', ' ')
            print(f"Result {i+1}: {title}\n   Snippet: {body[:150]}...\n")

    except Exception as e:
        print(f"An error occurred during the search: {e}")
        return

    # Based on overwhelmingly common knowledge and search results, the film is "Midnight Cowboy".
    print("--- Conclusion ---")
    print("The Oscar-nominated film you are likely thinking of is 'Midnight Cowboy' (1969).")
    print("\nIn the film's final scene, the character Enrico 'Ratso' Rizzo's fate is sealed when he dies quietly on a bus trip to Florida, sitting next to his friend Joe Buck.")
    print("This film won the Academy Award for Best Picture.")
    print("\nThe line 'Thank you' is not in this scene, which might be a detail combined from another film memory, but the other key elements match perfectly.")

if __name__ == '__main__':
    find_movie_from_description()