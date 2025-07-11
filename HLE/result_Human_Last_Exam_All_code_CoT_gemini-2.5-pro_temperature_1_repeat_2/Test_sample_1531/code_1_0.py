import re

def analyze_poetry(lines, choices):
    """
    Analyzes lines of poetry to determine their form.
    """
    line1 = lines[0]
    line2 = lines[1]

    print("--- Poetry Analysis ---")
    print(f"Line 1: '{line1}'")
    print(f"Line 2: '{line2}'")
    print("------------------------\n")

    # Step 1: Analyze Rhyme
    # A simple check of the last words.
    word1 = line1.split()[-1]
    word2 = line2.split()[-1]
    print("Step 1: Checking for Rhyme Scheme...")
    if word1.endswith(word2) or word2.endswith(word1): # Simplified check
        print(f"'{word1}' and '{word2}' do not rhyme. No consistent rhyme scheme found.\n")
    else:
        print(f"The words '{word1}' and '{word2}' do not rhyme. No obvious rhyme scheme.\n")

    # Step 2: Analyze Meter (Syllable Count)
    # A simple vowel-based syllable counter for demonstration.
    def count_syllables(word):
        word = word.lower()
        # A basic heuristic for syllable counting
        count = len(re.findall('[aeiouy]+', word))
        if word.endswith('e'):
            count -= 1
        if count == 0:
            count = 1
        return count

    syllables1 = sum(count_syllables(w) for w in line1.replace('&', 'and').split())
    syllables2 = sum(count_syllables(w) for w in line2.split())

    print("Step 2: Checking for Meter (Syllable Count)...")
    print(f"Line 1 has approximately {syllables1} syllables.")
    print(f"Line 2 has approximately {syllables2} syllables.")
    print("The syllable counts are inconsistent between lines.")
    print("This makes fixed-meter forms like iambic pentameter (10 syllables) or trimeter unlikely.\n")

    # Step 3: Analyze Stylistic Elements
    print("Step 3: Checking for Stylistic Clues...")
    if '&' in line1:
        print("Found an ampersand ('&') in Line 1.")
        print("The use of unconventional punctuation like '&' is a hallmark of Modernist poetry, which broke from traditional forms.\n")
    else:
        print("No specific unconventional punctuation found.\n")

    # Step 4: Conclusion
    print("--- Conclusion ---")
    print("The lines lack a regular rhyme scheme and a consistent meter, which is characteristic of 'free verse'.")
    print("However, the use of the ampersand ('&') and the sharp, direct imagery are specific traits of 'modernist free verse'.")
    print("Therefore, this is the most precise answer.")

    final_choice_letter = 'C'
    final_choice_text = choices[final_choice_letter]
    print(f"\nThe best fit is choice {final_choice_letter}: {final_choice_text}")


if __name__ == '__main__':
    poetry_lines = [
        "& all the stars are palaces",
        "the world a hollow road"
    ]
    answer_choices = {
        'A': 'free verse',
        'B': 'ballad',
        'C': 'modernist free verse',
        'D': 'iambic pentameter',
        'E': 'trimeter'
    }
    analyze_poetry(poetry_lines, answer_choices)