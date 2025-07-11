import requests
import re
import collections

def solve_chaucer_rhyme():
    """
    This script analyzes "The Book of the Duchess" to find which of a given
    list of words is not used in a rhyme.
    """
    print("Analyzing Chaucer's 'The Book of the Duchess' for rhyming words...\n")

    # 1. Fetch the text
    print("Step 1: Fetching the text from Project Gutenberg...")
    try:
        url = "https://www.gutenberg.org/cache/epub/2248/pg2248.txt"
        response = requests.get(url)
        response.raise_for_status()
        full_text = response.text
        print("Text fetched successfully.")
    except requests.exceptions.RequestException as e:
        print(f"Could not fetch text from URL: {e}. Please check your internet connection.")
        return

    # 2. Isolate the poem's text
    print("\nStep 2: Isolating the poem's text...")
    try:
        start_marker = "I have gret wonder, be this lyght,"
        end_marker = "Now goth goodly, for I wol, slepe."
        start_index = full_text.index(start_marker)
        end_index = full_text.index(end_marker) + len(end_marker)
        poem_text = full_text[start_index:end_index]
        print("Poem text isolated.")
    except ValueError:
        print("Could not find start/end markers. The text format may have changed.")
        return

    # 3. Process lines to get last words
    print("\nStep 3: Extracting the last word from each line...")
    lines = poem_text.splitlines()
    last_words = []
    for line in lines:
        # Find all word-like sequences in the line
        words = re.findall(r'\b\w+\b', line.lower())
        if words:
            last_words.append(words[-1])
    print(f"Extracted {len(last_words)} last words from the poem's lines.")

    # 4. Identify rhyming couplets
    print("\nStep 4: Identifying rhyming pairs from couplets...")
    rhyme_pairs = []
    # The poem is in AABB couplets, so we check pairs of consecutive lines
    for i in range(0, len(last_words) - 1, 2):
        word1 = last_words[i]
        word2 = last_words[i+1]
        rhyme_pairs.append((word1, word2))
    print(f"Identified {len(rhyme_pairs)} potential rhyming pairs.")

    # 5. Check target words
    print("\nStep 5: Checking each target word against the rhyming pairs...")
    target_words = ["wente", "here", "fool", "hool", "countour"]
    rhyme_findings = collections.defaultdict(list)
    word_in_text = {word: False for word in target_words}
    poem_text_lower = poem_text.lower()

    # First, a simple check if the word exists in the poem at all
    for word in target_words:
        if re.search(r'\b' + re.escape(word) + r'\b', poem_text_lower):
            word_in_text[word] = True

    # Now, check for rhymes
    for word1, word2 in rhyme_pairs:
        if word1 in target_words:
            # To be a rhyme, the partner word should be different
            if word1 != word2:
                rhyme_findings[word1].append(word2)
        if word2 in target_words:
            if word1 != word2:
                rhyme_findings[word2].append(word1)

    # 6. Report findings
    print("\n--- Analysis Results ---")
    not_found_in_rhyme = None
    for word in target_words:
        cap_word = word.capitalize()
        if not word_in_text[word]:
            print(f"Analysis for '{cap_word}': The word does not appear in the poem.")
            not_found_in_rhyme = cap_word
        elif rhyme_findings[word]:
            partners = sorted(list(set(rhyme_findings[word])))
            print(f"Analysis for '{cap_word}': Found rhyming with {partners}.")
        else:
            print(f"Analysis for '{cap_word}': Appears in the text but was not found in a rhyming couplet.")
            not_found_in_rhyme = cap_word

    print("\n--- Conclusion ---")
    if not_found_in_rhyme:
        print(f"The word that Chaucer does NOT make a rhyme with is '{not_found_in_rhyme}'.")
    else:
        print("All target words were found in rhyming pairs.")

if __name__ == '__main__':
    solve_chaucer_rhyme()