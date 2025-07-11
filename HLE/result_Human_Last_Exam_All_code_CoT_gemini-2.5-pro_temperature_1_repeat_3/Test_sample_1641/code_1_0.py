import collections

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3.
    """
    max_len = 3
    alphabet = ['a', 'b']

    # 1. Generate all words of length 0 to max_len
    words = []
    q = collections.deque([""])
    words.append("")
    while q:
        word = q.popleft()
        if len(word) < max_len:
            for char in alphabet:
                new_word = word + char
                words.append(new_word)
                q.append(new_word)
    
    # Although the prompt mentions 255 substitutions, the number of pairs (x,y)
    # with length <= 3 is 15 * 15 = 225. We will proceed with this number.
    
    # 2. Iterate through all pairs and count finite ones based on the rules.
    finite_len_lt = 0  # Counter for cases where len(y) < len(x)
    finite_len_eq = 0  # Counter for cases where len(y) = len(x)
    finite_len_gt = 0  # Counter for cases where len(y) > len(x)

    for x in words:
        # If x is the empty word, the substitution is non-finite.
        if not x:
            continue

        for y in words:
            # Case 1: len(y) < len(x) -> Finite
            if len(y) < len(x):
                finite_len_lt += 1
            
            # Case 2: len(y) == len(x) -> Finite if x != y
            elif len(y) == len(x):
                if x != y:
                    finite_len_eq += 1
            
            # Case 3: len(y) > len(x) -> Finite if x is not a subword of y
            elif len(y) > len(x):
                if x not in y:
                    finite_len_gt += 1

    total_finite = finite_len_lt + finite_len_eq + finite_len_gt

    # 3. Print the result as an equation.
    print(f"Number of finite substitutions:")
    print(f"{finite_len_lt} (for len(y) < len(x)) + {finite_len_eq} (for len(y) = len(x)) + {finite_len_gt} (for len(y) > len(x)) = {total_finite}")

solve()
<<<175>>>