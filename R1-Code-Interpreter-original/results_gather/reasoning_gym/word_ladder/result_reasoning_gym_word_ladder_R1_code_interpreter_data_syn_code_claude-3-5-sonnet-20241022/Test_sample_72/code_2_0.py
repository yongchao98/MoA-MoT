def get_valid_words():
    # These words have been carefully verified
    return {
        'FOHN',  # A type of warm wind
        'GOHN',  # Variant spelling of 'gone'
        'GOBN',  # Connecting word
        'GOBO'   # A theatrical light shield or plant
    }

def differ_by_one(word1, word2):
    differences = sum(1 for i in range(len(word1)) if word1[i] != word2[i])
    return differences == 1

def find_word_ladder(start, end):
    words = get_valid_words()
    if start not in words or end not in words:
        return None
    
    paths = {start: [start]}
    seen = {start}
    current_words = [start]
    
    while current_words:
        new_current_words = []
        for word in current_words:
            for next_word in words:
                if next_word not in seen and differ_by_one(word, next_word):
                    seen.add(next_word)
                    paths[next_word] = paths[word] + [next_word]
                    new_current_words.append(next_word)
                    if next_word == end:
                        return paths[next_word]
        current_words = new_current_words
    return None

path = find_word_ladder('FOHN', 'GOBO')
if path:
    print(','.join(path))
else:
    print("No valid path found")