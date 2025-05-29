def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word into a tuple of characters
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[sorted_word].append(word)
    
    # Return the list of anagram groups
    return list(anagrams.values())

words = ["eldern", "lender", "fp", "pf", "cohere", "cheero", "anglican", "canaling", "converters", "reconverts", "wandy", "dawny"]
result = group_anagrams(words)
print(result)