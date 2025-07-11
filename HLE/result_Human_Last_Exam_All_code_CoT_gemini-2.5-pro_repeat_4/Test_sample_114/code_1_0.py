import requests
import re

def solve_chaucer_rhyme_puzzle():
    """
    Analyzes Chaucer's "Book of the Duchess" to find which word from a given
    list is not used at the end of a line to make a rhyme.
    """
    # 1. Fetch the text from Project Gutenberg
    try:
        url = "https://www.gutenberg.org/cache/epub/2247/pg2247.txt"
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        full_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text: {e}")
        return

    # 2. Isolate the poem's text
    # The poem starts with 'I HAVE gret wonder...' and ends before the Gutenberg license.
    # A reliable end marker is the last line of the poem itself.
    start_phrase = "I HAVE gret wonder, be this lyght,"
    end_phrase = "Was hoolly in my remembraunce."
    
    try:
        start_index = full_text.index(start_phrase)
        end_index = full_text.index(end_phrase) + len(end_phrase)
        poem_text = full_text[start_index:end_index]
    except ValueError:
        print("Could not find the start or end of the poem in the downloaded text.")
        return

    # 3. Extract all words used at the end of lines
    lines = poem_text.splitlines()
    rhyming_words = set()
    
    for line in lines:
        # Find all words in the line
        words_in_line = re.findall(r'\b\w+\b', line)
        if words_in_line:
            # The last word is the one used for rhyming
            last_word = words_in_line[-1].lower()
            rhyming_words.add(last_word)

    # 4. Check the answer choices against the set of rhyming words
    choices = {
        'A': 'Wente',
        'B': 'Here',
        'C': 'Fool',
        'D': 'Hool',
        'E': 'Countour'
    }
    
    not_found_word = None
    not_found_letter = None

    print("Analyzing which word is not used in a rhyme in 'Book of the Duchess':\n")
    for letter, word in choices.items():
        word_lower = word.lower()
        if word_lower in rhyming_words:
            print(f"- The word '{word}' is used in a rhyme.")
        else:
            print(f"- The word '{word}' was NOT found as a rhyming word.")
            not_found_word = word
            not_found_letter = letter
    
    # 5. Output the final answer
    if not_found_word:
        print(f"\nConclusion: Chaucer does not make a rhyme with the word '{not_found_word}'.")
        print(f"<<<{not_found_letter}>>>")
    else:
        print("\nConclusion: All words were found in rhyming positions.")
        print("<<<F>>>")


if __name__ == "__main__":
    solve_chaucer_rhyme_puzzle()