import requests
import re

def solve_chaucer_rhyme_puzzle():
    """
    This function determines which word from a given list Chaucer did not make a
    rhyme with in "Book of the Duchess" by checking for the word's presence in the text.
    """
    # Define the target words and their corresponding answer choices.
    target_words = {
        'A': 'wente',
        'B': 'here',
        'C': 'fool',
        'D': 'hool',
        'E': 'countour'
    }

    # URL for a plain text version of "Book of the Duchess" from Project Gutenberg.
    url = "https://www.gutenberg.org/cache/epub/39531/pg39531.txt"

    print(f"Downloading text from {url}...")
    try:
        response = requests.get(url)
        # Raise an error if the download fails.
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text file. {e}")
        return

    # Isolate the main body of the poem to avoid metadata.
    try:
        start_marker = "I have gret wonder, be this light,"
        end_marker = "*** END OF THE PROJECT GUTENBERG EBOOK"
        start_index = text.index(start_marker)
        end_index = text.index(end_marker)
        poem_text = text[start_index:end_index]
    except ValueError:
        print("Error: Could not find the start/end markers of the poem in the downloaded file.")
        # As a fallback, use the whole text.
        poem_text = text

    # Clean the text: convert to lowercase and find all whole words.
    words_in_text = re.findall(r'\b[a-z]+\b', poem_text.lower())
    unique_words = set(words_in_text)

    print("\n--- Checking for each word in 'Book of the Duchess' ---")

    final_answer_key = 'F'  # Corresponds to "None of the above"
    not_found_word = None

    # Check for the presence of each target word.
    for key, word in sorted(target_words.items()):
        if word in unique_words:
            print(f"Result for '{word.capitalize()}' ({key}): FOUND")
        else:
            print(f"Result for '{word.capitalize()}' ({key}): NOT FOUND")
            final_answer_key = key
            not_found_word = word.capitalize()

    print("------------------------------------------------------")

    # Display the final conclusion.
    if not_found_word:
        print(f"\nConclusion: The word '{not_found_word}' was not found in the text.")
        print("Therefore, Chaucer did not make a rhyme with it in this poem.")
    else:
        print("\nConclusion: All words were found in the text.")

if __name__ == '__main__':
    solve_chaucer_rhyme_puzzle()