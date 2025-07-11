import sys
import io

# A simple syllable counting dictionary for the words in the poem
# This is a common method for this type of analysis, as programmatic
# syllabification can be complex.
syllable_counts = {
    'the': 1, 'first': 1, 'step': 1,
    'thirty': 2, 'five': 1,
    'rules': 1, 'and': 1, 'lines': 1,
    'an': 1, 'intricate': 3, "spider's": 2, 'web': 1, 'work': 1
}

# The full text of the poem, including title and page number
# We treat "thirty-five" as two words for counting purposes.
poem_words = [
    "The", "first", "step",
    "thirty", "five",
    "rules", "and", "lines",
    "an", "intricate", "spider's", "web", "work"
]

def analyze_poem_syllables():
    """
    Analyzes the syllable count of the poem and prints the breakdown.
    """
    total_syllables = 0
    equation_parts = []

    print("Analyzing the syllables of the poem:")
    print("-" * 35)

    # Process title words
    title_words = poem_words[:3]
    title_syllables = sum(syllable_counts[w.lower()] for w in title_words)
    print(f"Title ('The first step'): {' + '.join([str(syllable_counts[w.lower()]) for w in title_words])} = {title_syllables} syllables")
    total_syllables += title_syllables
    equation_parts.extend([str(syllable_counts[w.lower()]) for w in title_words])


    # Process page number words
    page_num_words = poem_words[3:5]
    page_num_syllables = sum(syllable_counts[w.lower()] for w in page_num_words)
    print(f"Page Number ('thirty five'): {' + '.join([str(syllable_counts[w.lower()]) for w in page_num_words])} = {page_num_syllables} syllables")
    total_syllables += page_num_syllables
    equation_parts.extend([str(syllable_counts[w.lower()]) for w in page_num_words])

    # Process collaged words
    collaged_words = poem_words[5:]
    collaged_syllables = sum(syllable_counts[w.lower()] for w in collaged_words)
    print(f"Collaged words: {' + '.join([str(syllable_counts[w.lower()]) for w in collaged_words])} = {collaged_syllables} syllables")
    total_syllables += collaged_syllables
    equation_parts.extend([str(syllable_counts[w.lower()]) for w in collaged_words])

    print("-" * 35)
    print("Final Calculation:")
    # Print the full equation for total syllables
    final_equation = " + ".join(equation_parts)
    print(f"Total syllables = {final_equation} = {total_syllables}")
    print("-" * 35)

    if total_syllables == 17:
        print("The poem has 17 syllables, which is the exact length of an American Sentence.")
    else:
        print(f"The poem has {total_syllables} syllables.")

analyze_poem_syllables()
