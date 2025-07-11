import itertools

def generate_words(alphabet, max_len):
    """
    Generates all words on the given alphabet up to a maximum length,
    including the empty word.
    """
    words = [""]
    for length in range(1, max_len + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
    return words

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite based on established criteria
    from string rewriting theory.
    """
    # Case 1: x is the empty word.
    # The substitution "" -> y can be applied infinitely if y is not empty.
    if not x:
        # It's finite only if y is also empty (the "" -> "" rule does nothing).
        return not y

    # Case 2: The length of y is less than or equal to the length of x.
    # The length of the string cannot grow, so the process must terminate.
    if len(y) <= len(x):
        return True

    # Case 3: The length of y is greater than the length of x.
    # The substitution is infinite if a new instance of x can be generated.
    # This happens if y contains x, or if x and y can overlap to recreate x.

    # a) y contains x as a subword.
    # e.g., x='a', y='bab'. Start with 'a', get 'bab', which contains 'a'. Infinite.
    if x in y:
        return False

    # b) A non-empty suffix of x is a prefix of y.
    # e.g., x='ab', y='baba'. Suffix 'b' of x is prefix 'b' of y.
    # A word like 'cab' becomes 'c(baba)', which contains 'ab'. Infinite.
    for i in range(1, len(x) + 1):
        suffix_x = x[-i:]
        if y.startswith(suffix_x):
            return False
            
    # c) A non-empty prefix of x is a suffix of y.
    # e.g., x='ab', y='cab'. Prefix 'ab' of x is suffix 'ab' of y.
    # A word like 'abab' becomes '(cab)ab', which contains 'ab'. Infinite.
    for i in range(1, len(x) + 1):
        prefix_x = x[:i]
        if y.endswith(prefix_x):
            return False

    # If none of the conditions for infinitude are met, the substitution is finite.
    return True

def solve_and_print():
    """
    Calculates the number of finite substitutions for words of length <= 3
    on the alphabet {a, b} and prints the detailed result.
    """
    alphabet = ['a', 'b']
    max_len = 3
    
    words = generate_words(alphabet, max_len)
    
    # The number of pairs (x,y) is 15*15=225.
    
    finite_counts_by_len_x = {0: 0, 1: 0, 2: 0, 3: 0}
    
    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_counts_by_len_x[len(x)] += 1
    
    total_finite = sum(finite_counts_by_len_x.values())
    
    print("Number of finite substitutions, broken down by length of x:")
    
    c0 = finite_counts_by_len_x[0]
    c1 = finite_counts_by_len_x[1]
    c2 = finite_counts_by_len_x[2]
    c3 = finite_counts_by_len_x[3]
    
    print(f"For x = '' (length 0): {c0}")
    print(f"For x of length 1: {c1}")
    print(f"For x of length 2: {c2}")
    print(f"For x of length 3: {c3}")
    print("\nFinal equation:")
    print(f"{c0} + {c1} + {c2} + {c3} = {total_finite}")
    print(f"\nTotal number of finite substitutions: {total_finite}")

solve_and_print()
<<<164>>>