def group_anagrams(words):
    anagrams = {}
    for word in words:
        # Sort the word to form the key
        key = ''.join(sorted(word))
        # Add the word to the corresponding anagram group
        if key in anagrams:
            anagrams[key].append(word)
        else:
            anagrams[key] = [word]
    # Return the grouped anagrams
    return list(anagrams.values())

words = ["unsober", "unrobes", "bournes", "pitsaw", "sawpit", "aristotle", "literatos", 
         "striolate", "thaness", "hastens", "snathes", "arid", "raid", "dari", 
         "dangering", "gandering", "gardening", "deranging"]

result = group_anagrams(words)
print(result)