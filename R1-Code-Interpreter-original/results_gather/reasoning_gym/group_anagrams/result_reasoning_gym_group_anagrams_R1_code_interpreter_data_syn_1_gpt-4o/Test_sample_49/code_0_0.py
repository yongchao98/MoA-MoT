def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word and use it as a key
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["sturnine", "runniest", "set", "est", "neurad", "endura", "unread", 
         "hotches", "shochet", "rehallow", "hallower", "dismalities", 
         "dissimilate", "halfungs", "flashgun", "tirelessness", 
         "ritelessness", "dirca", "carid", "daric", "darci", "batiste", "bastite"]

result = group_anagrams(words)
print(result)