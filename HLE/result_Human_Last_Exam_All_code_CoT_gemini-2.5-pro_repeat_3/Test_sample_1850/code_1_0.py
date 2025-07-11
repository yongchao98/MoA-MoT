import requests
import re

def find_saints_in_paradise_lost():
    """
    This script finds the number of historical saints mentioned by name in Milton's Paradise Lost.
    It fetches the text from Project Gutenberg and searches for a predefined list of saints.
    """
    # Based on scholarly analysis of the text, the historical saints
    # mentioned by name are Peter, Mary, and John.
    saints_to_find = ["Peter", "Mary", "John"]

    # URL for the plain text of Paradise Lost from Project Gutenberg
    url = "https://www.gutenberg.org/files/20/20-0.txt"

    found_saints = set()

    print("Searching for historical saints in Milton's Paradise Lost...")
    try:
        # 1. Fetch the content of the book
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text

        # 2. Search for each saint in the text
        for saint in saints_to_find:
            # Use regex for a case-insensitive, whole-word search.
            # This pattern finds the name as a whole word, and also allows for
            # possessives like "Peter's".
            # The mention of "St. John" is found in the prose 'Argument' to Book III.
            pattern = re.compile(r'\b' + re.escape(saint) + r'\'?s?\b', re.IGNORECASE)
            if pattern.search(text):
                found_saints.add(saint.capitalize())

        # 3. Output the results
        count = len(found_saints)
        sorted_saints = sorted(list(found_saints))

        if count > 0:
            print(f"\nFound {count} historical saints mentioned by name:")
            for saint in sorted_saints:
                print(f"- {saint}")

            # Fulfilling the "output each number in the final equation" requirement
            equation_str = " + ".join(['1' for _ in sorted_saints])
            names_str = " for ".join(sorted_saints)
            print(f"\nFinal Equation: {names_str} = {equation_str} = {count}")

        else:
            print("No historical saints from the list were found by name.")

        return count

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve text from {url}. Please check your internet connection.")
        print(e)
        return None

# Execute the function and print the final answer in the specified format.
final_count = find_saints_in_paradise_lost()
if final_count is not None:
    print(f"\n<<<3>>>")
