import itertools

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3.
    """
    # 1. Generate all words of length 0, 1, 2, 3 on the alphabet {a,b}
    alphabet = ['a', 'b']
    all_words = []
    for length in range(4): # Corresponds to lengths 0, 1, 2, 3
        # The empty word for length 0
        if length == 0:
            all_words.append("")
            continue
        # Generate all words of a given length
        for p in itertools.product(alphabet, repeat=length):
            all_words.append("".join(p))
    
    total_pairs = 0
    non_finite_count = 0

    # 2. Iterate through all possible couples (x, y)
    for x in all_words:
        for y in all_words:
            total_pairs += 1
            
            # 3. Check if the substitution x -> y is non-finite
            is_non_finite = False
            # Condition for non-finite: x is empty, or y is longer/same-length and contains x.
            if x == "":
                is_non_finite = True
            elif len(y) >= len(x):
                if x in y:
                    is_non_finite = True

            if is_non_finite:
                non_finite_count += 1

    # 4. Calculate the number of finite substitutions and print the result
    finite_count = total_pairs - non_finite_count
    
    print(f"Total possible substitutions: {total_pairs}")
    print(f"Number of non-finite substitutions: {non_finite_count}")
    print(f"Number of finite substitutions is the difference:")
    print(f"{total_pairs} - {non_finite_count} = {finite_count}")

solve()

# The final answer is the number of finite substitutions
print("<<<163>>>")