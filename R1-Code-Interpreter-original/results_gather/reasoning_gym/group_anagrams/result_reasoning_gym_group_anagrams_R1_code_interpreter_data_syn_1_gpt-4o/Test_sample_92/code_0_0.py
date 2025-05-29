from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to form the key
        key = ''.join(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[key].append(word)
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["hydronephrosis", "nephrohydrosis", "boarts", "tabors", "aborts", 
         "bonaci", "bianco", "alnico", "oilcan", "nonskier", "einkorns", 
         "cliftonia", "fictional"]

# Group the anagrams and print the result
result = group_anagrams(words)
print(result)