from collections import deque

# A simple dictionary of valid words for demonstration purposes
valid_words = set([
    "SHED", "SHAD", "SHID", "SHIB", "SHEB", "SHEA", "SHEE", "SHEF", "SHEG", "SHEH", "SHEI", "SHEJ", "SHEK", "SHEL", "SHEM", "SHEN", "SHEO", "SHEP", "SHEQ", "SHER", "SHES", "SHET", "SHEU", "SHEV", "SHEW", "SHEX", "SHEY", "SHEZ",
    "CHED", "CHID", "CHIB", "CRIB", "CRID", "CRAB", "CRAD", "CRIB", "CRIC", "CRID", "CRIE", "CRIF", "CRIG", "CRIH", "CRII", "CRIJ", "CRIK", "CRIL", "CRIM", "CRIN", "CRIO", "CRIP", "CRIQ", "CRIR", "CRIS", "CRIT", "CRIU", "CRIV", "CRIW", "CRIX", "CRIY", "CRIZ"
])

def is_valid_word(word):
    return word in valid_words

def word_ladder(start, end):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current_word[:i] + c + current_word[i+1:]
                if next_word not in visited and is_valid_word(next_word):
                    visited.add(next_word)
                    queue.append((next_word, path + [next_word]))
    
    return []

# Find the word ladder from SHED to CRIB
ladder = word_ladder("SHED", "CRIB")
print(",".join(ladder))