import collections

def find_number_anagram(s):
    """
    Checks if a string s is an anagram of a number word from "one" to "ten".
    Returns the number as an integer if an anagram is found, otherwise None.
    """
    s_counts = collections.Counter(s.lower())
    
    number_map = {
        "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
        "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10
    }
    
    for word, number in number_map.items():
        if len(s) == len(word):
            word_counts = collections.Counter(word)
            if s_counts == word_counts:
                return word, number
    return None, None

def decode_words():
    """
    Applies the decoding logic to the list of words and prints the result.
    """
    words = [
        "utpshtheas",
        "fkuhu",
        "numhkatnatu",
        "anhaeketn"
    ]
    
    print("Finding the pattern by checking the last three letters of each word...\n")
    
    final_number = None
    final_word_str = ""

    for word in words:
        last_three = word[-3:]
        number_word, number_val = find_number_anagram(last_three)
        
        if number_word:
            print(f"'{word}' -> last three letters are '{last_three}'.")
            print(f"'{last_three}' is an anagram for the word '{number_word}'.")
            print(f"The decoded number is {number_val}.\n")
            final_number = number_val
            final_word_str = number_word
        else:
            print(f"'{word}' -> last three letters are '{last_three}'. No number anagram found.\n")

    if final_number is not None:
        # The prompt requires printing each number in the final equation.
        # Here, the equation is the decoded number itself.
        print("The final decoded word is 'ten'.")
        print("The final decoded number is:")
        for digit in str(final_number):
            print(digit)

decode_words()
