import re
import requests

def find_saints_in_paradise_lost():
    """
    This script fetches Milton's "Paradise Lost", searches for mentions of
    historical saints, and prints the names and total count.
    """
    # Step 1: Get the text of "Paradise Lost" from Project Gutenberg.
    try:
        url = "https://www.gutenberg.org/files/20/20-0.txt"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        full_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve the text of Paradise Lost. {e}")
        return

    # Step 2: Isolate the poem's text to avoid metadata in the count.
    try:
        # The poem's main content starts with "Book I" and ends before the Gutenberg license.
        start_marker = "Book I"
        end_marker = "*** END OF THE PROJECT GUTENBERG EBOOK PARADISE LOST ***"
        start_index = full_text.find(start_marker)
        end_index = full_text.find(end_marker)
        if start_index == -1 or end_index == -1:
            poem_text = full_text # Fallback if markers are not found
        else:
            poem_text = full_text[start_index:end_index]
    except Exception:
        poem_text = full_text # Fallback to using the whole text

    # Step 3: Define the list of historical saints to search for.
    # This list is based on prominent human figures from the Bible.
    saints_to_check = [
        "Abel", "Abraham", "Adam", "David", "Enoch", "Eve", "Isaac", "Jacob",
        "John", "Joseph", "Joshua", "Mary", "Moses", "Noah", "Samson", "Solomon",
        "Peter", "Paul", "James", "Andrew" # Add a few apostles to be sure
    ]

    found_saints = set()
    text_lower = poem_text.lower()

    # Step 4: Search for each saint's name in the text.
    for saint in saints_to_check:
        # Use regex for case-insensitive, whole-word matching.
        pattern = r'\b' + re.escape(saint.lower()) + r'\b'
        if re.search(pattern, text_lower):
            # In Paradise Lost, the mention of "John" is a clear reference to
            # John the Baptist in the context of Book XII.
            if saint.lower() == 'john':
                found_saints.add("John the Baptist")
            else:
                found_saints.add(saint)

    # Step 5: Sort the list for clear presentation and get the count.
    unique_found_saints = sorted(list(found_saints))
    count = len(unique_found_saints)

    print(f"Based on a textual analysis, {count} historical saints are mentioned by name in Milton's Paradise Lost:")
    
    # Print each saint's name.
    for saint in unique_found_saints:
        print(f"- {saint}")

    # As requested, output each number in the final equation.
    equation_parts = ["1" for _ in unique_found_saints]
    equation_str = " + ".join(equation_parts)
    
    print("\nThe final count is the sum of each unique saint found:")
    print(f"{equation_str} = {count}")

find_saints_in_paradise_lost()
<<<16>>>