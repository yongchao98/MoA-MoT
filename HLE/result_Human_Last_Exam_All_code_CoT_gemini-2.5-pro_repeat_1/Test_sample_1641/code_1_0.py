import collections

def solve():
    """
    Calculates the number of finite substitutions x -> y for words x, y of length <= 3.
    """
    # 1. Generate all words of length 0 to 3
    words = [""]
    queue = [""]
    for _ in range(3):
        new_queue = []
        for word in queue:
            for char in ['a', 'b']:
                new_word = word + char
                words.append(new_word)
                new_queue.append(new_word)
        queue = new_queue

    # Group words by length for structured output
    words_by_len = collections.defaultdict(list)
    for word in words:
        words_by_len[len(word)].append(word)

    # 2. Iterate through all pairs and count finite substitutions
    # A substitution x -> y is finite iff x is not a substring of y.
    
    counts_by_len = []
    total_finite_count = 0
    
    # Iterate through lengths of x from 0 to 3
    for length in sorted(words_by_len.keys()):
        finite_count_for_len = 0
        for x in words_by_len[length]:
            for y in words:
                if x not in y:
                    finite_count_for_len += 1
        
        print(f"Number of finite substitutions for |x|={length}: {finite_count_for_len}")
        counts_by_len.append(str(finite_count_for_len))
        total_finite_count += finite_count_for_len

    # 3. Output the final equation and answer
    equation = " + ".join(counts_by_len)
    print(f"Total number of finite substitutions: {equation} = {total_finite_count}")

solve()
<<<164>>>