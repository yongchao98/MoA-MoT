import re

def solve_russian_stress():
    """
    Finds and prints Russian words from a given text that have more than one syllable
    and where the stress is not on the first syllable.
    """
    # The text to be analyzed
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary to store stress information. The key is the lowercase word,
    # and the value is the 0-indexed position of the stressed syllable.
    # For example, in "Шашлык" (Шаш-лы́к), the stress is on the 2nd syllable, so its index is 1.
    stress_data = {
        "шашлык": 1,   # Шаш-лы́к
        "запах": 0,    # за́-пах
        "горелым": 1,  # го-ре́-лым
        "вскоре": 0,   # вско́-ре
        "прибежал": 2, # при-бе-жа́л
        "сосед": 1,   # со-се́д
        "видно": 0,   # ви́д-но
        "дошёл": 1,    # до-шёл
        "его": 1,      # е-го́
        "дома": 0      # до́-ма (in the context of "to his house")
    }

    # Russian vowels, used for counting syllables
    vowels = "аеёиоуыэюя"

    # Split the text into words using spaces and hyphens as delimiters
    # and remove punctuation.
    words = re.split(r'[\s-]+', text)
    
    # List to store the final words that match the criteria
    result_words = []

    for word in words:
        # Clean the word of any non-alphabetic characters and convert to lowercase
        clean_word = re.sub(r'[^а-яА-ЯёЁ]', '', word).lower()

        if not clean_word:
            continue

        # Count syllables by counting the vowels
        syllable_count = sum(1 for char in clean_word if char in vowels)

        # Check if the word has more than one syllable and is in our stress dictionary
        if syllable_count > 1 and clean_word in stress_data:
            # Get the stressed syllable index
            stress_position = stress_data[clean_word]

            # If the stress is not on the first syllable (index > 0)
            if stress_position > 0:
                # Add the original word (with original casing) to our results
                result_words.append(re.sub(r'[^а-яА-ЯёЁ]', '', word))

    # Print the final list, comma-separated
    print(", ".join(result_words))

solve_russian_stress()