import requests
import re

def find_chaucer_rhymes():
    """
    Analyzes Chaucer's "Book of the Duchess" to find which word from a given
    list is not used in a rhyme.
    """
    # 1. Fetch the text from Project Gutenberg
    try:
        url = "https://www.gutenberg.org/cache/epub/2246/pg2246.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text: {e}")
        return

    # 2. Isolate the poem and extract the last word of each line
    start_marker = "I HAVE gret wonder, be this lyght,"
    end_marker = "Now goth goodly, for I wol, slepe."
    try:
        start_index = text.index(start_marker)
        # Find the final line and include it
        end_index = text.index(end_marker) + len(end_marker)
        poem_text = text[start_index:end_index]
    except ValueError:
        print("Could not find the start or end of the poem in the downloaded text.")
        return

    lines = poem_text.splitlines()
    last_words = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Split by spaces and get the last part
        words_in_line = line.split()
        if words_in_line:
            # Clean the word: lowercase, remove punctuation except internal hyphens
            last_word = re.sub(r'[^\w-]$', '', words_in_line[-1]).lower()
            last_words.append(last_word)

    choices = {
        "A": "wente",
        "B": "here",
        "C": "fool", # The poem uses the spelling 'fole'
        "D": "hool",
        "E": "countour"
    }

    rhyme_map = {}
    # Build a map of all rhyming pairs assuming couplets (word at even index rhymes with odd)
    for i in range(0, len(last_words) - 1, 2):
        word1 = last_words[i]
        word2 = last_words[i+1]
        rhyme_map[word1] = word2
        rhyme_map[word2] = word1
        
    # Handle spelling variation for "fool" -> "fole"
    if "fole" in rhyme_map:
        rhyme_map["fool"] = rhyme_map["fole"]


    print("Checking which word does not have a rhyme in 'Book of the Duchess':\n")
    unrhymed_choice = "F. None of the above"

    for key, word in choices.items():
        if word in rhyme_map:
            partner = rhyme_map[word]
            print(f"Choice {key} ('{word}') has a rhyme.")
            # Final output needs each number in the equation.
            # In this context, it's each word in the rhyme.
            print(f"Found rhyme: {word} / {partner}")
            print("-" * 20)
        else:
            # Check if the word appears in the text at all, but not at a line end
            if f" {word} " in poem_text.lower():
                 print(f"Choice {key} ('{word}') is NOT rhymed.")
                 print(f"The word '{word}' appears in the text but not at the end of a line, so it cannot be part of a rhyme.")
            else:
                 print(f"Choice {key} ('{word}') is NOT rhymed.")
                 print(f"The word '{word}' was not found at the end of any line.")

            unrhymed_choice = key
            # We don't need to break, show analysis for all
            print("-" * 20)

    print(f"\nThe word Chaucer does NOT make a rhyme with is '{choices.get(unrhymed_choice)}'.")


if __name__ == "__main__":
    find_chaucer_rhymes()
<<<E>>>