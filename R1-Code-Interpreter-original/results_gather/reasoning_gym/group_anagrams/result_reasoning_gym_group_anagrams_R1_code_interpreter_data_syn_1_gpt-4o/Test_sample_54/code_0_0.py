words = ["hor", "rho", "unremote", "routemen"]

# Dictionary to hold the groups of anagrams
anagrams = {}

# Group words by sorted tuple of characters
for word in words:
    # Sort the word to form the key
    key = tuple(sorted(word))
    # Add the word to the corresponding anagram group
    if key in anagrams:
        anagrams[key].append(word)
    else:
        anagrams[key] = [word]

# Extract the groups of anagrams
result = list(anagrams.values())

print(result)