import itertools

def solve_and_count():
    """
    Finds the number of finite substitutions for words of length <= 3 on the alphabet {a,b}.
    """
    
    # 1. Generate all words of length <= 3.
    words = []
    for length in range(4):  # Lengths 0, 1, 2, 3
        for p in itertools.product('ab', repeat=length):
            words.append("".join(p))
            
    # Total number of possible substitutions (x,y)
    total_substitutions = len(words) * len(words)
    non_finite_count = 0

    # 2. Iterate through all possible pairs (x, y).
    for x in words:
        for y in words:
            
            # 3. Check if the substitution x -> y is non-finite.
            is_non_finite = False
            
            # Rule for empty word x: non-finite if y is not empty.
            if not x:
                if y:
                    is_non_finite = True
            # Rule for |y| > |x|: check for "pumping" conditions.
            elif len(y) > len(x):
                # Condition 1: x is a subword of y.
                if x in y:
                    is_non_finite = True
                else:
                    # Condition 2: An overlap `s+y` contains `x`, where `s` is a proper suffix of `x`.
                    for i in range(1, len(x)):
                        suffix = x[i:]
                        if x in (suffix + y):
                            is_non_finite = True
                            break
                    
                    if not is_non_finite:
                        # Condition 3: An overlap `y+p` contains `x`, where `p` is a proper prefix of `x`.
                        for i in range(1, len(x)):
                            prefix = x[:i]
                            if x in (y + prefix):
                                is_non_finite = True
                                break
            
            # Increment count if found to be non-finite.
            if is_non_finite:
                non_finite_count += 1
    
    finite_count = total_substitutions - non_finite_count
    
    # 4. Print the numbers of the final equation as requested.
    print(total_substitutions, "-", non_finite_count, "=", finite_count)

solve_and_count()