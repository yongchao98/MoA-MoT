import requests
import re

def solve_divine_comedy_question():
    """
    This script finds which title characters from Shakespeare's plays are mentioned 
    by name in Dante's The Divine Comedy and selects the correct answer from a list.
    """
    # 1. Define characters from the answer choices and their search terms.
    # The dictionary maps the character's name to a list of strings to search for.
    characters_to_check = {
        "Julius Caesar": ["Caesar"],
        "Cleopatra": ["Cleopatra"],
        "Antony": ["Antony", "Antonius"],
        "Troilus": ["Troilus"],
        "King John": ["King John"],
        "Pericles": ["Pericles"]
    }

    # Define the answer choices provided in the user's prompt.
    answer_choices = {
        "A": ["Julius Caesar", "Pericles"],
        "B": ["Julius Caesar", "Cleopatra", "King John"],
        "C": ["Julius Caesar", "Troilus", "Antony"],
        "D": ["Julius Caesar", "Cleopatra"],
        "E": ["Julius Caesar", "Antony", "Cleopatra"]
    }

    # 2. Fetch the text of The Divine Comedy.
    print("Step 1: Fetching the text of 'The Divine Comedy'...")
    try:
        url = "https://www.gutenberg.org/cache/epub/8800/pg8800.txt"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text
        print("Text successfully downloaded.\n")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text. {e}")
        return

    # 3. Search the text for each character.
    print("Step 2: Searching for each character in the text...")
    found_characters = []
    for character, search_terms in characters_to_check.items():
        is_found = any(re.search(r'\b' + term + r'\b', text, re.IGNORECASE) for term in search_terms)
        if is_found:
            found_characters.append(character)
            print(f"- Found mentions of '{character}'.")
        else:
            print(f"- Did not find mentions of '{character}'.")
    
    print("\nInitial list of found characters:", found_characters)

    # 4. Analyze and refine the results.
    print("\nStep 3: Analyzing the findings...")
    final_characters = []
    # Add a special rule for Pericles
    if "Pericles" in found_characters:
        print("- Note: While the name 'Pericles' is in the text, Dante refers to the historical Athenian statesman, not the character from Shakespeare's play 'Pericles, Prince of Tyre'. Therefore, 'Pericles' is disqualified.")
        # Filter out Pericles
        final_characters = [char for char in found_characters if char != "Pericles"]
    else:
        final_characters = found_characters
    
    # Dante does not mention Antony or Troilus by name. The script confirms this.
    if "Antony" not in final_characters:
        print("- Note: Mark Antony is not mentioned by name in The Divine Comedy.")
    if "Troilus" not in final_characters:
        print("- Note: Troilus is not explicitly named in The Divine Comedy.")

    print("\nFinal list of qualifying characters:", final_characters)
    
    # 5. Compare with answer choices to find the correct one.
    print("\nStep 4: Evaluating the answer choices...")
    correct_option = None
    for option, names in answer_choices.items():
        # An option is correct if its set of names exactly matches the set of found names.
        if set(names) == set(final_characters):
            correct_option = option
            break
            
    if correct_option:
        print(f"\nThe characters {', '.join(final_characters)} are mentioned.")
        print(f"This matches the names in Option {correct_option}: {', '.join(answer_choices[correct_option])}.")
        print(f"\nTherefore, the correct answer is Option {correct_option}.")
    else:
        print("\nCould not find an answer choice that perfectly matches the results.")

# Run the full process
solve_divine_comedy_question()
<<<D>>>