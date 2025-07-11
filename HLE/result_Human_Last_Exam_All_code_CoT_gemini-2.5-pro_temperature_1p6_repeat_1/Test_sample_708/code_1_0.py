import re

def find_critical_word_position():
    """
    Analyzes a passage to find the word position where elevated reading times
    are expected due to a garden-path ambiguity.
    """
    passage = "Experience will be a decisive factor in this battle. The ships will outmaneuver the enemy fleet because the old man the boats."

    # Use regex to split the passage into words, handling punctuation.
    # \b matches a word boundary, \w+ matches one or more word characters.
    words = re.findall(r'\b\w+\b', passage)

    # The ambiguous phrase is "the old man". We are looking for the word "man"
    # when it acts as a verb following the subject "the old".
    target_word = "man"
    preceding_words = ["the", "old"]
    
    # We will find the ordinal position of the word "man" in the critical phrase.
    # We iterate through the words list to find the sequence "the", "old", "man".
    critical_word_position = -1
    for i in range(len(words) - 2):
        # Check if the current slice of words matches our target phrase.
        # We compare in lowercase to ensure the match is case-insensitive.
        if words[i].lower() == preceding_words[0] and \
           words[i+1].lower() == preceding_words[1] and \
           words[i+2].lower() == target_word:
            
            # The position is the index + 1 (since lists are 0-indexed).
            # We are interested in the position of "man", which is at index i+2.
            critical_word_position = i + 2 + 1
            break

    print("The passage is: \"{}\"".format(passage))
    print("\nThe ambiguity occurs in the phrase '...the old man the boats.'")
    print("The word 'man' is ambiguous. It can be a noun (part of 'the old man') or a verb (to staff).")
    print("Processing the less common verb interpretation causes elevated reading times, even when the reader is successful.")
    print("The critical word is therefore '{}'.".format(target_word))
    print("\nCalculating the ordinal word position of '{}':".format(target_word))
    
    # In our script, the calculation is finding the index in the word list.
    # The equation is: position = index_of_the_word + 1
    # Index of "man" is critical_word_position - 1.
    print("Found '{}' at index {} in the word list.".format(target_word, critical_word_position - 1))
    print("The final equation for ordinal position is: Position = Index + 1")
    print("Position = {} + 1 = {}".format(critical_word_position - 1, critical_word_position))
    print("\nThe final answer for the ordinal position is: {}".format(critical_word_position))

find_critical_word_position()