import requests
import re

def find_unrhymed_word():
    """
    Finds which word from a list does not appear in Chaucer's "Book of the Duchess".
    A word cannot be rhymed if it's not in the text.
    """
    # URL for the plain text of "Book of the Duchess" from Project Gutenberg
    url = "https://www.gutenberg.org/cache/epub/2247/pg2247.txt"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text file: {e}")
        return

    choices = {
        'A': "Wente",
        'B': "Here",
        'C': "Fool",
        'D': "Hool",
        'E': "Countour"
    }

    not_found_word_choice = None
    not_found_word = ""

    print("Checking for the presence of each word in 'The Book of the Duchess':\n")

    for choice_key, word in choices.items():
        # Use regex to find the whole word, case-insensitively.
        # \b ensures we match whole words only.
        pattern = r'\b' + re.escape(word) + r'\b'
        if re.search(pattern, text, re.IGNORECASE):
            print(f"Result for '{word}': Found in the text.")
        else:
            print(f"Result for '{word}': NOT found in the text.")
            if not_found_word_choice is None:
                not_found_word_choice = choice_key
                not_found_word = word
    
    if not_found_word_choice:
        print(f"\nThe word '{not_found_word}' is not in the text, so Chaucer could not have made a rhyme with it.")
    else:
        print("\nAll words were found in the text. The analysis might require phonetic checks.")

if __name__ == '__main__':
    find_unrhymed_word()
