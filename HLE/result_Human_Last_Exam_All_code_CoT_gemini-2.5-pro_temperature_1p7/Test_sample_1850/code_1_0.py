import requests
import re

def count_saints_in_paradise_lost():
    """
    Downloads Milton's "Paradise Lost" and counts the number of historical saints
    mentioned by name.
    """
    # URL for the plain text version of Paradise Lost from Project Gutenberg
    url = "https://www.gutenberg.org/files/20/20-0.txt"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text: {e}")
        return

    # List of historical saints (primarily biblical figures) to search for.
    # The Archangels (Michael, Gabriel, Raphael, Uriel) are considered saints in many traditions.
    saints_to_find = [
        "Adam", "Eve", "Enoch", "Noah", "Abraham", "Moses", "David", "Mary",
        "Michael", "Gabriel", "Raphael", "Uriel"
    ]
    
    found_saints = set()

    # Search for each saint in the text
    for saint in saints_to_find:
        # Use regex to find the name as a whole word, case-insensitively
        # \b ensures we match whole words only (e.g., 'Eve' not 'evening')
        pattern = r'\b' + re.escape(saint) + r'\b'
        if re.search(pattern, text, re.IGNORECASE):
            found_saints.add(saint)
    
    # Sort the found saints for consistent output
    sorted_found_saints = sorted(list(found_saints))
    count = len(sorted_found_saints)

    print("Historical saints mentioned by name in Milton's 'Paradise Lost':")
    if count > 0:
        # The user requested the "final equation", which we interpret as showing the components of the count.
        equation_str = " + ".join(["1" for _ in sorted_found_saints])
        names_str = ", ".join(sorted_found_saints)
        print(f"Found: {names_str}")
        print(f"Count: {equation_str} = {count}")
    else:
        print("None were found from the list.")

if __name__ == "__main__":
    count_saints_in_paradise_lost()