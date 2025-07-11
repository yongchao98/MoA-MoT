import itertools

def get_words_by_length():
    """Generates all words of length 0 to 3 on the alphabet {a, b}."""
    words_by_length = {i: [] for i in range(4)}
    words_by_length[0] = [""]
    for length in range(1, 4):
        for p in itertools.product(['a', 'b'], repeat=length):
            words_by_length[length].append("".join(p))
    return words_by_length

def solve():
    """
    Calculates the number of finite substitutions (x, y) for words of length <= 3.
    """
    words_by_length = get_words_by_length()
    all_words = sum(words_by_length.values(), [])
    
    # According to our analysis, a substitution x->y is finite if and only if
    # x is not the empty string, and one of the following conditions holds:
    # 1. len(y) < len(x)
    # 2. len(y) == len(x) AND x != y
    # 3. len(y) > len(x) AND x is not a substring of y

    # Category 1: len(y) < len(x)
    # x cannot be the empty string, as len(x) > len(y) >= 0.
    count1 = 0
    for len_x in range(1, 4):
        for len_y in range(len_x):
            count1 += len(words_by_length[len_x]) * len(words_by_length[len_y])

    # Category 2: len(y) == len(x) and x != y
    # x cannot be the empty string here.
    count2 = 0
    for length in range(1, 4):
        num_words = len(words_by_length[length])
        # Total pairs of this length is num_words * num_words
        # Pairs where x == y is num_words
        count2 += num_words * num_words - num_words

    # Category 3: len(y) > len(x) and x not in y
    # x cannot be the empty string.
    count3 = 0
    non_empty_words = all_words[1:]
    for x in non_empty_words:
        for y in all_words:
            if len(y) > len(x) and x not in y:
                count3 += 1

    total_finite = count1 + count2 + count3
    
    print(f"Finite substitutions where len(y) < len(x): {count1}")
    print(f"Finite substitutions where len(y) == len(x) and x != y: {count2}")
    print(f"Finite substitutions where len(y) > len(x) and x is not a subword of y: {count3}")
    print(f"Total finite substitutions = {count1} + {count2} + {count3} = {total_finite}")

solve()
<<<163>>>