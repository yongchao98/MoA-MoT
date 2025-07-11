import re
# The 'russian_stresses' library is required for this script to work.
# If not installed, you can install it by running: pip install russian_stresses
try:
    from russian_stresses import Accent
except ImportError:
    print("Error: The 'russian_stresses' library is not installed.")
    print("Please install it by running: pip install russian_stresses")
    exit()

def find_words_with_non_initial_stress(text):
    """
    Finds and prints words from a Russian text that meet two criteria:
    1. They have more than one syllable.
    2. The stress is not on the first syllable.
    The words are printed comma-separated and in their order of appearance.
    """
    # Define the set of Russian vowels for syllable counting.
    vowels = "аеёиоуыэюя"

    # Convert text to lowercase and extract all words using regex.
    # The pattern \b[а-яё]+\b ensures we get whole words containing Russian letters.
    words = re.findall(r'\b[а-яё]+\b', text.lower())

    # Initialize the accentuator to find word stresses.
    accentuator = Accent()

    result_words = []

    for word in words:
        # Count the number of vowels to determine the syllable count.
        syllable_count = sum(1 for char in word if char in vowels)

        # Condition 1: Exclude one-syllable words.
        if syllable_count <= 1:
            continue

        # Find the index of the first vowel in the word. This marks the first syllable.
        first_vowel_index = -1
        for i, char in enumerate(word):
            if char in vowels:
                first_vowel_index = i
                break

        stressed_vowel_index = -1
        # The letter 'ё' is always stressed in Russian.
        if 'ё' in word:
            stressed_vowel_index = word.find('ё')
        else:
            # Use the library to get the stressed version of the word.
            stressed_word = accentuator.put_accent(word)
            if stressed_word:
                # The stress mark '́' is added after the stressed vowel.
                stress_mark_index = stressed_word.find('́')
                if stress_mark_index != -1:
                    # The index of the stressed vowel in the original word
                    # is one less than the index of the stress mark.
                    stressed_vowel_index = stress_mark_index - 1
            else:
                # If the library cannot find the stress, we skip the word.
                continue

        # Condition 2: Check if the stress is NOT on the first syllable.
        # This is true if the stressed vowel's index is different from the first vowel's index.
        if stressed_vowel_index != -1 and stressed_vowel_index != first_vowel_index:
            result_words.append(word)

    # Print the final list of words, comma-separated.
    print(",".join(result_words))

# The input text provided by the user.
text_from_user = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

# Execute the function with the provided text.
find_words_with_non_initial_stress(text_from_user)