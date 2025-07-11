import requests
import re

def solve_saints_in_paradise_lost():
    """
    This function finds and counts the historical saints mentioned by name in
    Milton's Paradise Lost.
    """
    try:
        # 1. Fetch the text of Paradise Lost from Project Gutenberg
        url = "https://www.gutenberg.org/files/20/20-0.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text: {e}")
        return

    # 2. Isolate the text of Paradise Lost, including prose arguments
    try:
        # The work starts with the title.
        start_marker = "Paradise Lost."
        # It ends with a clear marker before Paradise Regained begins.
        end_marker = "THE END"
        start_index = text.find(start_marker)
        # Search for the end marker *after* the start marker
        end_index = text.find(end_marker, start_index)
        
        if start_index == -1 or end_index == -1:
            raise ValueError("Could not find start/end markers for Paradise Lost.")
            
        content = text[start_index:end_index]
    except ValueError as e:
        print(e)
        # Fallback to using the whole text if markers fail, though less accurate
        content = text

    # 3. Define a list of saints to search for
    saints_to_check = [
        "Peter", "Paul", "Mary", "Dominic", "Francis", "Augustine", "Jerome",
        "John", "Andrew", "James", "Philip", "Thomas", "Matthew", "Simon",
        "Jude", "Stephen", "Gregory", "Benedict"
    ]

    found_saints = set()

    # 4. Search for each saint's name in the text
    for saint in saints_to_check:
        # Use a regular expression for a case-insensitive, whole-word search.
        # For "Francis", use a negative lookahead to exclude "Franciscan".
        if saint.lower() == 'francis':
            pattern = re.compile(r'\bFrancis\b(? ![Cc][Aa][Nn]) ', re.IGNORECASE)
        else:
            pattern = re.compile(r'\b' + re.escape(saint) + r'\b', re.IGNORECASE)

        if pattern.search(content):
            found_saints.add(saint.capitalize())

    # 5. Print the results in the specified equation format
    if not found_saints:
        print("No historical saints were found by name in Paradise Lost.")
        print("0 = 0")
    else:
        # Sort the names for consistent output
        sorted_saints = sorted(list(found_saints))
        total = len(sorted_saints)
        
        # Create the equation string like "1 (Name1) + 1 (Name2)"
        equation_parts = [f"1 ({name})" for name in sorted_saints]
        equation_str = " + ".join(equation_parts)

        print(f"There are {total} historical saints mentioned by name in Milton's Paradise Lost.")
        print("The final count is derived as follows:")
        print(f"{equation_str} = {total}")

solve_saints_in_paradise_lost()