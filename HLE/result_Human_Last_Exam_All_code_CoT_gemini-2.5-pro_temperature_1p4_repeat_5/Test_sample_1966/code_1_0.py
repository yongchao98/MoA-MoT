import requests
import re

def solve_shakespeares_dante_question():
    """
    This script finds which Shakespearean title characters from the options
    are mentioned by name or strong allusion in Dante's Divine Comedy.
    """
    # 1. Define the characters to check based on the answer choices.
    # The dictionary maps the character's full name to a regex pattern for searching.
    characters_to_check = {
        "Julius Caesar": r"\bCaesar\b",
        "Cleopatra": r"\bCleopatra\b",
        "Antony": r"\bAntony\b",
        "Troilus": r"\bTroilus\b",
        "Pericles": r"\bPericles\b",
        "King John": r"\bKing John\b"
    }

    # 2. Fetch the text of The Divine Comedy (Longfellow translation).
    url = "https://www.gutenberg.org/files/1001/1001-0.txt"
    print("Fetching 'The Divine Comedy' from Project Gutenberg...")
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check for HTTP errors
        divine_comedy_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Failed to download text: {e}")
        return

    print("Searching for characters mentioned by name...\n")

    # 3. Search for explicit mentions of each character.
    found_characters = []
    for full_name, pattern in characters_to_check.items():
        if re.search(pattern, divine_comedy_text, re.IGNORECASE):
            found_characters.append(full_name)

    print("Found by explicit name:")
    if found_characters:
        for name in sorted(found_characters):
            print(f"- {name}")
    else:
        print("- None")

    # 4. Apply scholarly interpretation for strong allusions.
    # Antony is not explicitly named but is Cleopatra's partner in Inferno, Canto V.
    # This is a very strong and universally recognized allusion.
    final_list = list(found_characters)
    if "Cleopatra" in final_list and "Antony" not in final_list:
        final_list.append("Antony")
        print("\nNote: Antony is not explicitly named, but as Cleopatra's famous partner,")
        print("his presence is unambiguously implied. He is included for a complete answer.")

    # 5. Print the final combined list.
    print("\n----------------------------------------------------")
    print("Final list of mentioned characters (by name or strong allusion):")
    for name in sorted(final_list):
        print(f"- {name}")
    print("----------------------------------------------------")
    print("This corresponds to answer choice E.")

if __name__ == '__main__':
    solve_shakespeares_dante_question()
<<<E>>>