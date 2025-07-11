import urllib.request
import re
import string

def solve_chaucer_rhyme():
    """
    Analyzes Chaucer's "The Book of the Duchess" to find which
    of a given list of words is not used in a rhyming couplet.
    """
    # Step 1: Define target words from the multiple-choice question.
    target_words = {
        'A': 'wente',
        'B': 'here',
        'C': 'fool',
        'D': 'hool',
        'E': 'countour'
    }
    print("Task: Find the word Chaucer does NOT make a rhyme with from the following choices:")
    for key, word in target_words.items():
        print(f"  {key}. {word.capitalize()}")
    print("-" * 70)

    # Step 2: Download the text from Project Gutenberg.
    try:
        url = "https://www.gutenberg.org/files/2418/2418-0.txt"
        with urllib.request.urlopen(url) as response:
            # Use utf-8-sig to handle potential Byte Order Mark (BOM) at the file start
            raw_text = response.read().decode('utf-8-sig')
    except Exception as e:
        print(f"Error: Failed to download the text of the poem. Details: {e}")
        return

    # Step 3: Process the text to isolate the poem and get the last word of each line.
    try:
        start_marker = "I have so many an ydel thoght"
        end_marker = "End of the Project Gutenberg EBook"
        poem_text = raw_text.split(start_marker, 1)[1]
        poem_text = poem_text.split(end_marker, 1)[0]
    except IndexError:
        print("Error: Could not find the start/end markers to isolate the poem.")
        return

    all_lines = poem_text.strip().split('\r\n')
    line_data = []
    # Create a translator to remove all punctuation and digits. Hyphens are included for words like 'ful-cool'.
    translator = str.maketrans('', '', string.punctuation + string.digits)
    
    for line in all_lines:
        stripped_line = line.strip()
        if stripped_line:
            # Clean the line for word extraction: lowercase, no punctuation.
            clean_line = stripped_line.lower().translate(translator)
            words = clean_line.split()
            if words:
                line_data.append({'original': stripped_line, 'last_word': words[-1]})

    # Step 4: Define a simple function to check for rhymes.
    def is_rhyme(w1, w2):
        """A simple heuristic for Middle English rhyme: checks for common endings."""
        if not w1 or not w2:
            return False
        # For this poem, checking if the last 3 letters match is a reasonable heuristic.
        return w1[-3:] == w2[-3:]

    # Step 5 & 6: Analyze each target word to see if it's used in a rhyme.
    print("Analyzing each target word against the poem...\n")
    final_answer_key = None

    for key, word_to_check in target_words.items():
        print(f"Checking word: '{word_to_check.capitalize()}'")
        is_rhymed_in_poem = False
        found_at_all = False
        
        # Search the entire poem for instances of the word at the end of a line.
        for i, data in enumerate(line_data):
            if data['last_word'] == word_to_check:
                found_at_all = True
                
                # Check the next line to see if it forms a rhyming couplet.
                if i < len(line_data) - 1:
                    next_data = line_data[i+1]
                    if is_rhyme(data['last_word'], next_data['last_word']):
                        is_rhymed_in_poem = True
                        print("  - Verdict: RHYME FOUND")
                        print(f"    Line: \"{data['original']}\"")
                        print(f"    Line: \"{next_data['original']}\"")
                        break # A rhyme was found; no need to check other occurrences.
        
        if not is_rhymed_in_poem and found_at_all:
            # The word was found, but never in a rhyming pair. This is our answer.
            final_answer_key = key
            print("  - Verdict: NO RHYME FOUND")
            # Show the non-rhyming instance from the text.
            for i, data in enumerate(line_data):
                if data['last_word'] == word_to_check:
                    print("    Example of a non-rhyming occurrence:")
                    print(f"    Line: \"{data['original']}\"")
                    if i < len(line_data) - 1:
                         print(f"    Line: \"{line_data[i+1]['original']}\"")
                    break

        print("-" * 30)

    # Step 7: Announce the final result based on the analysis.
    print("\n--- Final Conclusion ---")
    if final_answer_key:
        unrhymed_word = target_words[final_answer_key].capitalize()
        print(f"The analysis of the Project Gutenberg text shows that '{unrhymed_word}' appears")
        print("at the end of a line but does not form a recognizable rhyme with its partner line.")
        print(f"Therefore, the word Chaucer does NOT make a rhyme with is '{unrhymed_word}'.")
    else:
        print("Analysis suggests all words from the list are used in rhyming couplets in the poem.")

# Execute the main function to solve the problem.
solve_chaucer_rhyme()
<<<E>>>