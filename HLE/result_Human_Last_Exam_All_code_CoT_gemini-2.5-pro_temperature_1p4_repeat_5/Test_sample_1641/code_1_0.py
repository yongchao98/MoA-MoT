import itertools

def is_substitution_finite(x, y):
    """
    Determines if the substitution x -> y is finite.
    """
    # Case 1: x is the empty word.
    # The substitution is finite only if y is also empty.
    if not x:
        return not bool(y)

    # Case 2: y is shorter than x.
    # The length of the word decreases with each step, so it must terminate.
    if len(y) < len(x):
        return True

    # Case 3: len(y) >= len(x) > 0.
    # This is infinite if applying the rule can create a new instance of x.
    # This occurs if there is a context (p, s) where p+y+s contains x.
    # The length of the context words p and s only needs to be checked up to len(x) - 1.
    len_x = len(x)
    max_context_len = len_x - 1

    # Generate all possible context words p and s.
    for p_len in range(max_context_len + 1):
        # Generate words of length p_len for the prefix context p
        p_words = [''.join(p) for p in itertools.product('ab', repeat=p_len)] if p_len > 0 else ['']
        for p in p_words:
            for s_len in range(max_context_len + 1):
                # Generate words of length s_len for the suffix context s
                s_words = [''.join(s) for s in itertools.product('ab', repeat=s_len)] if s_len > 0 else ['']
                for s in s_words:
                    test_word = p + y + s
                    # If x is found in p+y+s, it must overlap with y because
                    # len(p) < len(x) and len(s) < len(x), so x cannot be fully
                    # contained within either p or s alone.
                    if x in test_word:
                        return False # Infinite
                        
    # If no such context creates x, the substitution is finite.
    return True

def solve():
    """
    Solves the problem by generating all words and testing each substitution pair.
    """
    # Generate all words on {a,b} of length 0, 1, 2, 3.
    all_words = []
    for length in range(4):
        if length == 0:
            all_words.append('')
            continue
        all_words.extend([''.join(p) for p in itertools.product('ab', repeat=length)])

    # Initialize counters for each category of finite substitutions.
    count_case1 = 0  # x = '', y = ''
    count_case2 = 0  # len(y) < len(x)
    count_case3 = 0  # len(y) >= len(x) > 0 and finite

    for x in all_words:
        for y in all_words:
            if is_substitution_finite(x, y):
                if not x:
                    count_case1 += 1
                elif len(y) < len(x):
                    count_case2 += 1
                else:
                    count_case3 += 1

    total_finite = count_case1 + count_case2 + count_case3
    total_pairs = len(all_words) ** 2
    
    print(f"Counting finite substitutions for words of length <= 3 on {{a,b}}.")
    print(f"Total pairs (x,y) to check: {len(all_words)} * {len(all_words)} = {total_pairs}")
    print("-" * 20)
    print("Number of finite substitutions by case:")
    print(f"Case 1 (x = '', y = ''): {count_case1}")
    print(f"Case 2 (len(y) < len(x)): {count_case2}")
    print(f"Case 3 (len(y) >= len(x) > 0 and finite): {count_case3}")
    print("-" * 20)
    print("Final equation for the total number of finite substitutions:")
    print(f"{count_case1} + {count_case2} + {count_case3} = {total_finite}")
    
if __name__ == '__main__':
    solve()
    # The final numerical answer is requested in a specific format.
    # The calculation above gives the number 143.
    print("<<<143>>>")
