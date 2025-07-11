import collections

def get_words():
    """Generates all words of length 0 to 3 on alphabet {a,b}."""
    words = [""]
    for i in range(3):
        words.extend([w + s for w in words if len(w) == i for s in ['a', 'b']])
    return words

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite based on established criteria.
    Returns a string code for the finiteness reason, or None if infinite.
    """
    if x == "":
        return "case_empty" if y == "" else None

    if len(y) < len(x):
        return "case_len_decrease"

    if len(y) == len(x):
        return "case_len_equal" if x != y else None

    # len(y) > len(x)
    if x in y:
        return None

    # Check for conjugates
    # A substitution x -> y is infinite if x = uv (u,v non-empty)
    # and y contains vu.
    if len(x) >= 2:
        for i in range(1, len(x)):
            u = x[:i]
            v = x[i:]
            conjugate = v + u
            if conjugate in y:
                return None
    
    return "case_len_increase_finite"

def solve():
    """
    Counts the number of finite substitutions for words of length <= 3.
    """
    words = get_words()
    
    # Check for prompt's number of substitutions, though it's likely a typo
    # num_substitutions = len(words) * len(words) # 15*15=225

    finite_counts = collections.defaultdict(int)

    for x in words:
        for y in words:
            finiteness_case = is_finite(x, y)
            if finiteness_case:
                finite_counts[finiteness_case] += 1
    
    # Deconstruct the counts for the final print statement
    c_empty = finite_counts["case_empty"]
    c_len_decrease = finite_counts["case_len_decrease"]
    c_len_equal = finite_counts["case_len_equal"]
    c_len_increase = finite_counts["case_len_increase_finite"]
    total_finite = c_empty + c_len_decrease + c_len_equal + c_len_increase

    print(f"Number of finite substitutions by case:")
    print(f"- Finiteness from x=\"\", y=\"\": {c_empty}")
    print(f"- Finiteness from |y|<|x|: {c_len_decrease}")
    print(f"- Finiteness from |y|=|x| and x!=y: {c_len_equal}")
    print(f"- Finiteness from |y|>|x| (special cases): {c_len_increase}")
    print(f"Total: {c_empty} + {c_len_decrease} + {c_len_equal} + {c_len_increase} = {total_finite}")

solve()

<<<159>>>