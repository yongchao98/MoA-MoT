from collections import deque

# A small set of valid words for demonstration purposes
valid_words = set([
    "DEAL", "DEIL", "DELL", "DEAL", "DEAR", "DEER", "DEAL", "DEAL", "DEAL", "DEAL",
    "WEAL", "WEAR", "WEIR", "WELL", "WILL", "WISP", "WELD", "WILD", "WISP", "WISP",
    "SEAL", "SEAR", "SEER", "SELL", "SILL", "SILK", "SIMP", "SIP", "SIR", "SIRP",
    "WISP"
])

def is_valid_word(word):
    return word in valid_words

def word_ladder(start, end):
    queue = deque([(start, [start])])
    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current_word[:i] + c + current_word[i+1:]
                if next_word != current_word and is_valid_word(next_word):
                    queue.append((next_word, path + [next_word]))
                    valid_words.remove(next_word)  # Avoid revisiting
    return []

# Find the word ladder
ladder = word_ladder("DEAL", "WISP")
print(",".join(ladder))