import itertools

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3
    on the alphabet {a, b}.
    """
    alphabet = ['a', 'b']
    max_len = 3
    
    # 1. Generate all words of length <= 3, including the empty word.
    words = [""]
    for length in range(1, max_len + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
            
    # Total number of pairs is 15 * 15 = 225.
    
    # 2. Initialize counters for each category of finite substitutions.
    # We categorize based on the conditions for finiteness.
    count_y_is_empty = 0
    count_len_y_less_than_x = 0
    count_len_y_equals_x = 0
    count_len_y_greater_than_x = 0

    # 3. Iterate through all pairs (x, y) and check for finiteness.
    for x in words:
        for y in words:
            len_x = len(x)
            len_y = len(y)

            # If x is the empty word, the substitution is non-finite.
            if len_x == 0:
                continue

            # Case 1: y is the empty word (and x is not). Finite.
            if len_y == 0:
                count_y_is_empty += 1
                continue

            # Case 2: len(y) < len(x). Finite.
            if len_y < len_x:
                count_len_y_less_than_x += 1
                continue

            # Case 3: len(y) == len(x). Finite if x != y.
            if len_y == len_x:
                if x != y:
                    count_len_y_equals_x += 1
                continue

            # Case 4: len(y) > len(x). Check for non-finite conditions.
            if len_y > len_x:
                is_finite = True
                # Condition a: x is a subword of y
                if x in y:
                    is_finite = False
                else:
                    # Condition b: Overlap (prefix of y is proper suffix of x)
                    for i in range(1, len(y) + 1):
                        prefix_y = y[:i]
                        if len(prefix_y) < len_x and x.endswith(prefix_y):
                            is_finite = False
                            break
                    if is_finite:
                        # Condition b: Overlap (suffix of y is proper prefix of x)
                        for i in range(len(y)):
                            suffix_y = y[i:]
                            if len(suffix_y) < len_x and x.startswith(suffix_y):
                                is_finite = False
                                break
                
                if is_finite:
                    count_len_y_greater_than_x += 1

    # 4. Calculate the total and print the breakdown.
    total_finite = count_y_is_empty + count_len_y_less_than_x + count_len_y_equals_x + count_len_y_greater_than_x
    
    print("Number of finite substitutions, broken down by case:")
    print(f"y is empty: {count_y_is_empty}")
    print(f"len(y) < len(x): {count_len_y_less_than_x}")
    print(f"len(y) == len(x) and x != y: {count_len_y_equals_x}")
    print(f"len(y) > len(x) and no growth conditions met: {count_len_y_greater_than_x}")
    print("---")
    print(f"Total: {count_y_is_empty} + {count_len_y_less_than_x} + {count_len_y_equals_x} + {count_len_y_greater_than_x} = {total_finite}")

solve()