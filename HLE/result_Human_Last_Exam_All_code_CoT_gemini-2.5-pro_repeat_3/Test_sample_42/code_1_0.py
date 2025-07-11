import math

def get_char_pos(char):
    """Returns the 1-based alphabetical position of a character."""
    return ord(char.lower()) - ord('a') + 1

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def find_word_that_breaks_pattern(words):
    """
    Finds the word in a list that breaks the pattern.
    Pattern: pos(L1) + pos(L3) - pos(L2) is NOT a perfect square.
    The function finds the word where the result IS a perfect square.
    """
    vowels = "aeiou"
    for word in words:
        if len(word) == 3 and word[0] not in vowels and word[1] in vowels and word[2] not in vowels:
            p1 = get_char_pos(word[0])
            p2 = get_char_pos(word[1])
            p3 = get_char_pos(word[2])
            
            result = p1 + p3 - p2
            
            if is_perfect_square(result):
                print(f"The word '{word}' does not follow the pattern.")
                print("The calculation is:")
                # The final print statement outputs each number in the equation.
                print(f"{p1} + {p3} - {p2} = {result}")
                print(f"The result, {result}, is a perfect square.")
                return word
    return None

# The list of words from the answer choices
answer_choices = ['leg', 'dam', 'rat', 'car', 'bin']

# Find and print the result
find_word_that_breaks_pattern(answer_choices)
