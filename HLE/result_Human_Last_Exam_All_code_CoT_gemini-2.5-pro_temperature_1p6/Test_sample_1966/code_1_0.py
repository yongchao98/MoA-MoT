import requests
import re

def solve_dante_shakespeare_question():
    """
    This script determines which Shakespearean title characters are mentioned
    in Dante's Divine Comedy by searching the original Italian text.
    """
    try:
        # 1. Fetch the source text (Dante's Divine Comedy in Italian)
        url = "https://www.gutenberg.org/files/1001/1001-0.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        response.encoding = 'utf-8'
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the source text. {e}")
        print("Cannot solve the problem without the text.")
        return

    # 2. Define characters and their Italian names to search for.
    # Note: Antony (Marco Antonio), Troilus (Troilo), Pericles (Pericle),
    # and King John (Giovanni Senzaterra) are not found in the text,
    # but we include them for a complete search.
    characters_to_search = {
        "Julius Caesar": r"\bCesare\b",
        "Cleopatra": r"\bCleopatràs\b", # Dante uses this specific form in Inferno V
        "Antony": r"\bMarco Antonio\b",
        "Troilus": r"\bTroilo\b",
        "Pericles": r"\bPericle\b",
        "King John": r"\bGiovanni Senzaterra\b"
    }

    print("Searching for Shakespearean title characters in Dante's *The Divine Comedy*...")
    print("-" * 30)

    found_characters = []
    # 3. Search the text for each character
    for name, pattern in characters_to_search.items():
        # Using re.IGNORECASE for case-insensitivity
        if re.search(pattern, text, re.IGNORECASE):
            found_characters.append(name)
            print(f"[FOUND]     {name}")
        else:
            print(f"[NOT FOUND] {name}")

    # 4. Report findings and identify the correct answer
    print("-" * 30)
    print("\nConclusion:")
    if not found_characters:
        print("No matching characters were found in the text.")
        return

    print("The characters from the choices mentioned by name in The Divine Comedy are:")
    for char in found_characters:
        print(f"- {char}")

    # 5. Determine which option matches the found characters
    options = {
        "A": ["Julius Caesar", "Pericles"],
        "B": ["Julius Caesar", "Cleopatra", "King John"],
        "C": ["Julius Caesar", "Troilus", "Antony"],
        "D": ["Julius Caesar", "Cleopatra"],
        "E": ["Julius Caesar", "Antony", "Cleopatra"],
    }

    correct_option = None
    for option_letter, option_list in options.items():
        if set(found_characters) == set(option_list):
            correct_option = option_letter
            break

    if correct_option:
        print(f"\nThis corresponds to answer choice {correct_option}.")
    else:
        print("\nNone of the answer choices perfectly match the findings.")

# Run the solver
solve_dante_shakespeare_question()
<<<D>>>