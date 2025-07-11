import itertools

def solve_word_substitutions():
    """
    This function calculates the number of finite substitutions x -> y
    where x and y are words on the alphabet {a, b} with length at most 3.
    It follows the plan of counting the infinite substitutions and subtracting from the total.
    """
    alphabet = ['a', 'b']
    max_len = 3

    # 1. Generate all words of length <= max_len.
    words = [""]
    for length in range(1, max_len + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
    
    # 2. Calculate the total number of possible substitutions.
    # The set contains 1 (len 0) + 2 (len 1) + 4 (len 2) + 8 (len 3) = 15 words.
    # Total pairs (x,y) are 15 * 15 = 225.
    total_substitutions = len(words) * len(words)
    
    # 3. Count the number of infinite substitutions.
    # A substitution x -> y is infinite if len(y) > len(x) AND y contains x.
    infinite_count = 0
    for x in words:
        for y in words:
            if len(y) > len(x) and x in y:
                infinite_count += 1
                
    # 4. The number of finite substitutions is the remainder.
    finite_count = total_substitutions - infinite_count
    
    # 5. Print the final calculation as an "equation".
    print(f"{total_substitutions} - {infinite_count} = {finite_count}")

solve_word_substitutions()
<<<177>>>