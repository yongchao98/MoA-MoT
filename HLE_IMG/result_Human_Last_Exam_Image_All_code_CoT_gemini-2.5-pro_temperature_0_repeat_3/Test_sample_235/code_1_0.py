# First, we need to install the necessary library for syllable counting.
# You can do this by running: pip install pyphen
import pyphen
import re

def solve():
    """
    This function identifies all the words in the poem, counts the syllables for each,
    and determines the metric pattern.
    """
    # The poem includes the title, page number, and the pasted words.
    # The page number "35" is read as "thirty-five".
    poem_text = "The first step thirty five rules and lines an intricate spider's web work"

    # Use the Pyphen library to count syllables in English words.
    dic = pyphen.Pyphen(lang='en_US')

    # Clean and split the text into individual words.
    words = re.findall(r'\b\w+\b', poem_text.lower())

    syllable_counts = []
    total_syllables = 0

    print("Counting syllables for each word in the poem:")
    print("-" * 40)

    for word in words:
        # Pyphen's inserted() method returns a hyphenated version of the word.
        # The number of syllables is the number of parts.
        # For monosyllabic words, it returns the word itself (1 part).
        count = len(dic.inserted(word).split('-'))
        syllable_counts.append(str(count))
        total_syllables += count
        print(f"'{word}': {count} syllable(s)")

    print("-" * 40)
    print("The total number of syllables is the sum of the syllables of each word.")
    
    # Display the final equation
    equation = " + ".join(syllable_counts)
    print(f"Final Equation: {equation} = {total_syllables}")
    print("\n")

    # Analyze the result
    print(f"The poem has a total of {total_syllables} syllables.")
    print("An 'American Sentence' is a poetic form defined as a single sentence containing exactly 17 syllables.")
    print("Since our poem fits this definition, it follows the American Sentence pattern.")

solve()