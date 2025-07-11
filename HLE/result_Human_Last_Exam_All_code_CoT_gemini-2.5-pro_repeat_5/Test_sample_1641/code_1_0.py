import collections

def is_finite(x, y):
    """
    Checks if the substitution x -> y is finite based on established criteria.

    A substitution is infinite if:
    1. x is the empty string and y is not.
    2. len(y) > len(x) and x is a substring of y.
    3. len(y) == len(x) and y is a cyclic permutation of x.
    
    All other substitutions are finite.
    """
    # Case for an empty 'x': infinite if 'y' is not empty (allows infinite insertions).
    if x == "":
        return y == ""

    # Length-decreasing substitutions are always finite.
    if len(y) < len(x):
        return True
    
    # Length-preserving substitutions: infinite if y is a cyclic permutation of x.
    # This check (x in y+y) cleverly handles all cyclic permutations, including x == y.
    elif len(y) == len(x):
        is_cyclic = x in (y + y)
        return not is_cyclic
        
    # Length-increasing substitutions: infinite if x is a substring of y.
    else:  # len(y) > len(x)
        return x not in y

def generate_words(alphabet, max_len):
    """
    Generates all unique words on the alphabet up to max_len, including the empty word.
    """
    words = {""}
    queue = collections.deque([""])
    while queue:
        current_word = queue.popleft()
        if len(current_word) < max_len:
            for char in alphabet:
                new_word = current_word + char
                words.add(new_word)
                queue.append(new_word)
    return list(words)

def solve():
    """
    Generates all substitutions and counts how many are finite.
    """
    alphabet = ['a', 'b']
    max_len = 3
    
    # 1. Generate all possible words for x and y.
    words = generate_words(alphabet, max_len)
    
    finite_count = 0
    total_count = 0
    
    # 2. Iterate through all possible pairs of (x, y) and check for finiteness.
    for x in words:
        for y in words:
            total_count += 1
            if is_finite(x, y):
                finite_count += 1
    
    infinite_count = total_count - finite_count
    
    # 3. Print the results.
    print(f"Total number of words (length <= {max_len}): {len(words)}")
    print(f"Total substitutions (x, y) checked: {total_count}")
    print(f"Number of infinite substitutions: {infinite_count}")
    print(f"Number of finite substitutions: {finite_count}")
    print("---")
    print(f"Final Equation: {total_count} (Total) - {infinite_count} (Infinite) = {finite_count} (Finite)")

solve()
<<<153>>>