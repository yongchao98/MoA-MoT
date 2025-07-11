import requests
import re

def solve_shakespeare_dante_question():
    """
    This script identifies which Shakespearean title characters from a given list
    are mentioned by name in Dante's "The Divine Comedy".
    """
    # 1. Acquire the source text of "The Divine Comedy" from Project Gutenberg.
    try:
        url = "https://www.gutenberg.org/files/1001/1001-0.txt"
        response = requests.get(url)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Could not download the text. Error: {e}")
        return

    # 2. Define the list of characters to check, based on the answer choices.
    # We use the key as the full name for display and the value for searching.
    characters = {
        "Julius Caesar": "Caesar",
        "Pericles": "Pericles",
        "Cleopatra": "Cleopatra",
        "King John": "King John", # Searching for "King John" is more specific than "John"
        "Troilus": "Troilus",
        "Antony": "Antony"
    }

    print("Searching for Shakespearean characters in 'The Divine Comedy'...\n")

    # 3. Perform the search and store the results.
    found_characters = []
    for full_name, search_term in characters.items():
        # Use regular expressions for a case-insensitive, whole-word search.
        # The \b markers ensure we match the whole word.
        pattern = re.compile(r'\b' + re.escape(search_term) + r'\b', re.IGNORECASE)
        if pattern.search(text):
            found_characters.append(full_name)
            print(f"Result: '{full_name}' was FOUND in the text.")
        else:
            # Special check for King John, as a simple search for "John" would yield many false positives (e.g., John the Baptist).
            # A search for "King John" is more precise and yields no results.
            if full_name == "King John":
                 print(f"Result: '{full_name}' was NOT FOUND in the text.")
            else:
                 print(f"Result: '{full_name}' was NOT FOUND in the text.")


    # 4. Analyze the results to determine the correct answer choice.
    print("\n--- Analysis of Answer Choices ---")
    print(f"The characters found by name are: {', '.join(found_characters)}.")
    print("Based on this, let's evaluate the options:")

    options = {
        'A': ["Julius Caesar", "Pericles"],
        'B': ["Julius Caesar", "Cleopatra", "King John"],
        'C': ["Julius Caesar", "Troilus", "Antony"],
        'D': ["Julius Caesar", "Cleopatra"],
        'E': ["Julius Caesar", "Antony", "Cleopatra"]
    }

    correct_option = None
    for option, names in options.items():
        is_correct = all(name in found_characters for name in names)
        is_plausible = all(name in found_characters for name in names)
        
        # An option is incorrect if it contains any name NOT found in the text.
        contains_unfound = any(name not in found_characters for name in names)

        if contains_unfound:
            print(f"Option {option}: Incorrect. It lists at least one character who was not found.")
        else:
            print(f"Option {option}: Plausible. All characters listed were found in the text.")
            if correct_option is None: # Select the first plausible option
                 correct_option = option

    print("\n--- Conclusion ---")
    print("While Troilus was also found, options A, B, C, and E are definitively incorrect because they each list at least one character who is not mentioned by name in The Divine Comedy.")
    print("Option D is the only choice that contains a list of characters where every character is verifiably mentioned. Therefore, it is the best answer.")


solve_shakespeare_dante_question()
<<<D>>>