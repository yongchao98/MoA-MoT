import re
from russian_stresses import Stresses

def find_words_with_non_initial_stress(text):
    """
    Finds words in a Russian text that are multi-syllable and not stressed on the first syllable.

    Args:
        text (str): The input Russian text.

    Returns:
        None: Prints the result directly.
    """
    # Instantiate the Stresses class. This may take a moment on first run
    # as it loads a dictionary of word stresses.
    try:
        stresses_db = Stresses()
    except Exception as e:
        print(f"Error initializing the russian-stresses library: {e}")
        print("Please make sure the library is installed with 'pip install russian-stresses'")
        return

    # A string containing all Russian vowels
    vowels = "аеёиоуыэюя"

    # Use regex to find all word sequences, converting to lowercase.
    # This handles punctuation and extracts clean words.
    words = re.findall(r'\b[а-яА-ЯёЁ-]+\b', text.lower())

    result_words = []
    processed_words = set()

    # Iterate through the words in their order of appearance
    for word in words:
        # Count the number of vowels (syllables) in the word
        num_vowels = sum(1 for char in word if char in vowels)

        # We are only interested in words with more than one syllable
        if num_vowels > 1:
            # Find the character index of the first vowel in the word
            first_vowel_index = -1
            for i, char in enumerate(word):
                if char in vowels:
                    first_vowel_index = i
                    break
            
            # This should always find a vowel because we already checked num_vowels > 1
            if first_vowel_index == -1:
                continue

            try:
                # Get stress information for the word from the library
                # The result is a list of possible stress patterns; we use the first (most common) one.
                stress_info = stresses_db.get_stresses(word, 'all')
                
                # Check if the library found any stress patterns for the word
                if stress_info:
                    # The index of the character that is the stressed vowel
                    stressed_vowel_char_index = stress_info[0]['stressed_vowel_index']

                    # Check if the stress is NOT on the first syllable
                    if stressed_vowel_char_index != first_vowel_index:
                        result_words.append(word)

            except Exception:
                # This block will catch errors if a word is not in the stress dictionary.
                # We can choose to ignore such words.
                pass

    # Print the final list of words, comma-separated.
    print(", ".join(result_words))

# The input text from the user
russian_text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

# Execute the function
find_words_with_non_initial_stress(russian_text)