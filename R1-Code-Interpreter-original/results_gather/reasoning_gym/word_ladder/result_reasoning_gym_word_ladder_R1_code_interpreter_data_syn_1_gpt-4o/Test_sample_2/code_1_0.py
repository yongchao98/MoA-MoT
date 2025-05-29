def is_valid_word(word, word_list):
    return word in word_list

def transform_word_ladder(start, end, word_list):
    from collections import deque

    def neighbors(word):
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != word[i]:
                    new_word = word[:i] + c + word[i+1:]
                    if is_valid_word(new_word, word_list):
                        yield new_word

    queue = deque([(start, [start])])
    visited = set([start])

    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        for neighbor in neighbors(current_word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None

# Define a small set of valid words
valid_words = {"TURD", "TURN", "BURN", "BARN", "BARD", "CARD", "CART", "BART", "BABA"}

# Find the transformation path
path = transform_word_ladder("TURD", "BABA", valid_words)
if path:
    print(",".join(path))
else:
    print("No valid transformation found.")