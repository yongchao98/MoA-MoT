def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word into a tuple of characters
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[sorted_word].append(word)
    
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

# List of words to group
words = ["tanproof", "antproof", "converser", "conserver", "sneaksman", "masskanne", "outlive", "ovulite", "stylate", "stately"]

# Group the anagrams and print the result
print(group_anagrams(words))