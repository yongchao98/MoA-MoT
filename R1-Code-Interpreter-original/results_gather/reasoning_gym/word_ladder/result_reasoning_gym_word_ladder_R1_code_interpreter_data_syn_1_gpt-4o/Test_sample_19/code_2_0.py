# Define a small set of words that might be useful for the transformation
word_list = {
    "SHIV", "SNAP", "SLAP", "TRAP", "TARP"
}

def is_valid_word(word):
    return word in word_list

def find_word_ladder(start, end):
    from collections import deque

    def neighbors(word):
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != word[i]:
                    new_word = word[:i] + c + word[i+1:]
                    if is_valid_word(new_word):
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

ladder = find_word_ladder('SHIV', 'TARP')
if ladder:
    print(','.join(ladder))
else:
    print("No valid word ladder found.")