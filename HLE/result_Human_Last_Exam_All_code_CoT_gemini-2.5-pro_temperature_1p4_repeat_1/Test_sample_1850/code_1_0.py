import requests
import re

def solve_milton_saints_query():
    """
    Finds and counts the number of historical saints mentioned by name
    in Milton's Paradise Lost.
    """
    # URL for the plain text of Paradise Lost from Project Gutenberg
    url = "https://www.gutenberg.org/files/20/20-0.txt"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        text = response.text

        # We are only interested in Milton's text, not the Project Gutenberg header/footer.
        start_marker = "Book I"
        end_marker = "*** END OF THE PROJECT GUTENBERG EBOOK PARADISE LOST ***"
        start_index = text.find(start_marker)
        end_index = text.find(end_marker)

        if start_index == -1 or end_index == -1:
            print("Could not find the start or end markers in the text.")
            return

        content = text[start_index:end_index]

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve the text of Paradise Lost. {e}")
        return

    # List of historical figures from the Old Testament generally considered saints.
    # Milton's prophetic vision given to Adam in Books XI and XII names these figures.
    # New Testament saints like Peter or Paul are not mentioned by name.
    saints_to_check = [
        "Adam", "Eve", "Abel", "Enoch", "Noah", "Abraham",
        "Isaac", "Jacob", "Joseph", "Moses", "Joshua", "David"
    ]

    found_saints = []
    # Use a set to store unique names found, then sort for consistent output
    unique_saints_found = set()

    for saint in saints_to_check:
        # Search for the name as a whole word (\b is a word boundary)
        # re.IGNORECASE makes the search case-insensitive
        if re.search(r'\b' + re.escape(saint) + r'\b', content, re.IGNORECASE):
            unique_saints_found.add(saint)

    found_saints = sorted(list(unique_saints_found))
    count = len(found_saints)

    if count > 0:
        print("The following historical saints are mentioned by name in Paradise Lost:")
        print(", ".join(found_saints))

        # As requested, output each number in the final equation
        equation_components = ["1"] * count
        equation_str = " + ".join(equation_components)

        print("\nThe final equation is:")
        print(f"{equation_str} = {count}")
    else:
        print("No historical saints from the list were found by name.")
        print("0")

solve_milton_saints_query()
<<<12>>>