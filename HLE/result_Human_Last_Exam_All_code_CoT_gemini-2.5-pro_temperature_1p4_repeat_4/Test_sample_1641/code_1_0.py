import itertools

def generate_words(max_len):
    """Generates all words on the alphabet {a,b} up to a given maximum length."""
    alphabet = ['a', 'b']
    words = [''] # Start with the empty word for length 0
    for length in range(1, max_len + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
    return words

def count_finite_substitutions():
    """
    Counts the number of finite substitutions for words of length <= 3.
    """
    words = generate_words(3)
    total_pairs = len(words) * len(words)
    finite_count = 0

    # These lists will store the pairs for the explanation in the output
    finite_pairs = []
    infinite_pairs = []

    for x in words:
        for y in words:
            is_finite = False
            # Apply the rules for finiteness
            if x == '':
                if y == '':
                    is_finite = True
            else: # x is not the empty string
                if len(y) < len(x):
                    is_finite = True
                elif len(y) == len(x):
                    if x != y:
                        is_finite = True
                else: # len(y) > len(x)
                    if x not in y:
                        is_finite = True
            
            if is_finite:
                finite_count += 1
                finite_pairs.append((x, y))
            else:
                infinite_pairs.append((x, y))

    print(f"Total number of possible words for x and y (length <= 3): {len(words)}")
    print(f"Total number of substitutions (x,y): {total_pairs}\n")
    print(f"Number of finite substitutions: {finite_count}")
    print(f"Number of infinite substitutions: {len(infinite_pairs)}")
    
    # Detailed breakdown for clarity, matching the manual calculation
    # Count infinite pairs as a cross-check
    inf_x_empty = sum(1 for y in words if y != '')
    
    inf_len_eq = sum(1 for w in words if w != '')

    inf_len_gt = 0
    for x in words:
        if x == '': continue
        for y in words:
            if len(y) > len(x) and x in y:
                inf_len_gt += 1
    
    manual_infinite_total = inf_x_empty + inf_len_eq + inf_len_gt
    
    print("\n--- Breakdown of Infinite Substitutions ---")
    print(f"1. x is empty, y is not empty: {inf_x_empty}")
    print(f"2. |y| = |x| and x = y (and x is not empty): {inf_len_eq}")
    print(f"3. |y| > |x| and y contains x: {inf_len_gt}")
    print(f"Total calculated infinite: {inf_x_empty} + {inf_len_eq} + {inf_len_gt} = {manual_infinite_total}")

    print(f"\n--- Final Answer ---")
    print(f"Out of {total_pairs} substitutions, the number of finite ones is: {finite_count}")


if __name__ == '__main__':
    count_finite_substitutions()
