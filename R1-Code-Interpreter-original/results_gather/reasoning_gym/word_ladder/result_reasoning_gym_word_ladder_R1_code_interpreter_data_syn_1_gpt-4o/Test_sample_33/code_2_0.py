def word_ladder(start, end, word_list):
    from collections import deque

    def is_valid_word(word):
        return word in word_list

    def get_neighbors(word):
        neighbors = []
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != word[i]:
                    new_word = word[:i] + c + word[i+1:]
                    if is_valid_word(new_word):
                        neighbors.append(new_word)
        return neighbors

    queue = deque([(start, [start])])
    visited = set([start])

    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path

        for neighbor in get_neighbors(current_word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return []

# A small dictionary of valid words for this problem
word_list = {"JIGS", "JEST", "JETS", "JITS", "JOTS", "JAGS", "JAMS", "JARS", "JAYS", "JIBS", "JIGS", "JINS", "JOBS", "JOGS", "JOSS", "JUST"}

# Find the word ladder
ladder = word_ladder("JIGS", "JEST", word_list)
print(",".join(ladder))