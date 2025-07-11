import requests
import re

def find_unrhymed_word():
    """
    Finds which word from a list Chaucer does not rhyme in "Book of the Duchess"
    by checking which words appear at the end of a line.
    """
    # 1. Download the text of "The Book of the Duchess".
    url = "https://www.gutenberg.org/files/2247/2247-0.txt"
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Failed to retrieve the text. An error occurred: {e}")
        return

    # 2. Isolate the poem's content.
    try:
        start_marker = "THE BOOK OF THE DUCHESS."
        end_marker = "*** END OF THE PROJECT GUTENBERG EBOOK THE BOOK OF THE DUCHESS ***"
        start_index = text.find(start_marker)
        end_index = text.find(end_marker)
        if start_index == -1 or end_index == -1:
            raise ValueError("Poem markers not found")
        poem_text = text[start_index:end_index]
    except ValueError as e:
        print(f"Could not isolate the poem from the text: {e}")
        # Fallback to using the whole text if markers are not found
        poem_text = text

    # 3. & 4. Process each line to create a set of cleaned end-words.
    lines = poem_text.splitlines()
    end_words = set()

    for line in lines:
        line_stripped = line.strip()
        if not line_stripped:
            continue
        
        words_in_line = line_stripped.split()
        if words_in_line:
            last_word = words_in_line[-1]
            # Clean the word: remove non-alphabetic characters and make it lowercase.
            cleaned_word = re.sub(r'[^a-zA-Z]', '', last_word).lower()
            if cleaned_word:
                end_words.add(cleaned_word)

    # 5. Check each target word against the set of end-words.
    target_words = {
        "A": "wente",
        "B": "here",
        "C": "fool",
        "D": "hool",
        "E": "countour"
    }

    unrhymed_word_details = None

    print("--- Analysis Report ---")
    print("Checking which words from the list appear at the end of a line in the poem.")
    print("-" * 20)
    
    for option, word in target_words.items():
        if word in end_words:
            print(f"Found: '{word.capitalize()}' (Option {option}) is used at the end of a line.")
        else:
            print(f"NOT Found: '{word.capitalize()}' (Option {option}) is NEVER used at the end of a line.")
            unrhymed_word_details = (option, word)

    print("-" * 20)

    # 6. Announce the result.
    print("\n--- Conclusion ---")
    if unrhymed_word_details:
        option, word = unrhymed_word_details
        print(f"The word that Chaucer does NOT make a rhyme with is '{word}'.")
        print(f"This is because the word '{word}' never appears in a rhyming position (at the end of a line) in the entire poem.")
        # This is commented out because the final answer should be in <<<>>> format below
        # print(f"The correct option is {option}.")
    else:
        print("All words from the list were found at the end of a line, which is unexpected.")
        print("The answer is likely F: None of the above.")

if __name__ == '__main__':
    find_unrhymed_word()