import collections

def get_words():
    """Generates all words of length 0 to 3 on alphabet {a, b}."""
    words = [""]
    for length in range(1, 4):
        current_len_words = []
        if length == 1:
            current_len_words = ["a", "b"]
        else:
            for w in words:
                if len(w) == length - 1:
                    current_len_words.append(w + "a")
                    current_len_words.append(w + "b")
        words.extend(current_len_words)
    return words

def are_conjugates(x, y):
    """Checks if two non-empty strings x and y are conjugates."""
    if len(x) != len(y) or len(x) == 0:
        return False
    # This is an efficient way to check for conjugacy
    return (x + x).find(y) != -1 and x != y

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3.
    """
    words = get_words()
    
    # The prompt states there are 255 substitutions, but the specified
    # constraints (words of length <= 3 on {a,b}) lead to 15*15=225 pairs.
    # We proceed based on the explicit constraints.
    total_pairs = len(words) * len(words)

    finite_counts_by_len = collections.defaultdict(int)
    
    # x is the word to be substituted
    for x in words:
        # If x is the empty string, the substitution is never finite.
        if x == "":
            continue

        # y is the word x is substituted with
        for y in words:
            # A substitution x -> y is finite if:
            # 1. x is not a substring of y
            # 2. x and y are not conjugates
            
            # Check for substring condition
            is_substring = x in y
            
            # Check for conjugacy condition
            # Note: The problem of non-termination for x=y is covered by the is_substring check.
            # We check for x != y in are_conjugates to avoid double counting x=y as a conjugate of itself.
            are_conj = are_conjugates(x, y)

            if not is_substring and not are_conj:
                finite_counts_by_len[len(x)] += 1
                
    total_finite_count = sum(finite_counts_by_len.values())
    
    print("Calculating the number of finite substitutions (x -> y):")
    len1_count = finite_counts_by_len[1]
    len2_count = finite_counts_by_len[2]
    len3_count = finite_counts_by_len[3]

    print(f"Number of finite pairs for |x|=1: {len1_count}")
    print(f"Number of finite pairs for |x|=2: {len2_count}")
    print(f"Number of finite pairs for |x|=3: {len3_count}")
    print(f"Total number of finite substitutions: {len1_count} + {len2_count} + {len3_count} = {total_finite_count}")

solve()
<<<148>>>