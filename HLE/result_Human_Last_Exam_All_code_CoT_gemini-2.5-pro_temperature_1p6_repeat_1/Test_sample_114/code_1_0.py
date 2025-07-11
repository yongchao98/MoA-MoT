import requests
import re

def solve_chaucer_rhyme_puzzle():
    """
    Analyzes Chaucer's "Book of the Duchess" to find which word
    from a given list is not present and therefore not rhymed.
    """
    # URL for the plain text of "The Book of the Duchess"
    url = "https://www.gutenberg.org/cache/epub/2245/pg2245.txt"
    words_to_check = ['wente', 'here', 'fool', 'hool', 'countour']
    word_options = {'wente': 'A', 'here': 'B', 'fool': 'C', 'hool': 'D', 'countour': 'E'}

    print("Analyzing 'The Book of the Duchess' by Geoffrey Chaucer...")
    print("----------------------------------------------------------")

    try:
        # 1. Fetch the text
        response = requests.get(url)
        response.raise_for_status()
        text = response.text.lower() # Convert to lowercase for easy searching

        # 2. Isolate the poem's content
        # The poem starts after the initial metadata. A good starting line is "I have gret wonder".
        # It ends before the Gutenberg license footer.
        start_marker = "i have gret wonder, be this light"
        end_marker = "*** end of the project gutenberg ebook the book of the duchess ***"
        
        start_index = text.find(start_marker)
        end_index = text.find(end_marker)

        if start_index == -1 or end_index == -1:
            print("Error: Could not find the start or end of the poem in the text.")
            return

        poem_text = text[start_index:end_index]
        
        # Use a regular expression to find whole words
        not_found_word = None
        not_found_option = None

        # 3. Check for the presence of each word
        print("Checking for the presence of each word in the text:\n")
        for word in words_to_check:
            # We search for the word as a whole word to avoid matching substrings
            # e.g., finding "here" in "there"
            if re.search(r'\b' + re.escape(word) + r'\b', poem_text):
                print(f"Word '{word.capitalize()}' (Option {word_options[word]}): Found in the text.")
            else:
                print(f"Word '{word.capitalize()}' (Option {word_options[word]}): NOT found in the text.")
                not_found_word = word.capitalize()
                not_found_option = word_options[word]
        
        print("----------------------------------------------------------")
        if not_found_word:
            print(f"Conclusion: The word Chaucer does NOT make a rhyme with is '{not_found_word}',")
            print("because it does not appear in 'The Book of the Duchess'.")
        else:
            print("Conclusion: All words were found in the text.")


    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching the text: {e}")

solve_chaucer_rhyme_puzzle()
<<<E>>>