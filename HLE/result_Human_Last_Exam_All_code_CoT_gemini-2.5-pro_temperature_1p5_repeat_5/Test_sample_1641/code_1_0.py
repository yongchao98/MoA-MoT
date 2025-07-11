import collections

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite based on established criteria.

    Args:
        x (str): The word to be replaced.
        y (str): The word to substitute in.

    Returns:
        bool: True if the substitution is finite, False otherwise.
    """
    len_x = len(x)
    len_y = len(y)

    # Condition 1: If x is the empty word, the substitution is infinite.
    if len_x == 0:
        return False

    # Condition 2: If the substitution is length-decreasing, it's always finite.
    if len_y < len_x:
        return True

    # Condition 3: If length-preserving, it's finite iff x is not y.
    if len_y == len_x:
        return x != y

    # Condition 4: If length-increasing, check for subwords and overlaps.
    if len_y > len_x:
        # 4a: If x is a subword of y, it's trivially infinite.
        if x in y:
            return False

        # 4b: Check if y starts with a proper suffix of x.
        # This can lead to an infinite sequence by "pumping" from the left.
        for i in range(1, len_x):
            suffix = x[i:]
            if y.startswith(suffix):
                return False

        # 4c: Check if y ends with a proper prefix of x.
        # This can lead to an infinite sequence by "pumping" from the right.
        for i in range(1, len_x):
            prefix = x[:i]
            if y.endswith(prefix):
                return False
        
        # If none of the infinite conditions are met, it's finite.
        return True
    
    return False # Should not be reached

def solve_and_print():
    """
    Calculates the number of finite substitutions for all possible
    couples (x, y) of words of length up to 3 on the alphabet {a, b}.
    """
    alphabet = ['a', 'b']
    max_len = 3
    
    # Step 1: Generate all words of length 0 to max_len
    words = ['']
    if max_len >= 1:
        q = collections.deque(list(alphabet))
        words.extend(q)
        current_len_words = list(q)
        for _ in range(max_len - 1):
            next_len_words = []
            for w in current_len_words:
                for char in alphabet:
                    next_word = w + char
                    next_len_words.append(next_word)
            words.extend(next_len_words)
            current_len_words = next_len_words

    # Step 2: Iterate through all pairs and count finite ones
    finite_count = 0
    total_pairs = len(words) * len(words)
    
    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_count += 1
    
    print(f"Total number of substitutions considered (for words of length <= 3): {total_pairs}")
    print(f"Number of finite substitutions found: {finite_count}")

# Execute the main function to find the answer
solve_and_print()
<<<128>>>