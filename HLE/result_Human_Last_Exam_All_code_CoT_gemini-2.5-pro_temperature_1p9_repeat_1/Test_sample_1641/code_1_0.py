import itertools

def count_finite_substitutions():
    """
    Calculates the number of finite substitutions for words of length <= 3
    on the alphabet {a, b}.
    """
    alphabet = ['a', 'b']
    words = ['']
    for length in range(1, 4):
        # Create all words of a given length
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
            
    # Total number of possible words of length 0, 1, 2, or 3.
    # 1 (len 0) + 2 (len 1) + 4 (len 2) + 8 (len 3) = 15 words.
    num_words = len(words)
    
    # Total possible substitutions (x, y) is num_words * num_words.
    # This equals 15 * 15 = 225. We assume 255 in the prompt is a typo.
    total_substitutions = num_words * num_words

    non_finite_count = 0
    
    # A substitution x -> y is non-finite if and only if:
    # 1. x is the empty string "" and y is not.
    # 2. x is a non-empty proper substring of y (i.e., x is in y and x != y).
    
    for x in words:
        for y in words:
            is_non_finite = False
            
            # Condition 1: x is empty, y is not.
            if len(x) == 0 and len(y) > 0:
                is_non_finite = True
            
            # Condition 2: x is a proper substring of y.
            # This implies len(x) < len(y). We also check len(x) > 0 to
            # avoid double-counting with condition 1.
            elif len(x) > 0 and x in y and x != y:
                is_non_finite = True

            if is_non_finite:
                non_finite_count += 1
                
    finite_count = total_substitutions - non_finite_count
    
    # The final output prints each number in the final equation.
    print(f"{total_substitutions} - {non_finite_count} = {finite_count}")

count_finite_substitutions()
<<<181>>>