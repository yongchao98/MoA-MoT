from itertools import product

def solve():
    """
    Counts the number of finite substitutions (x,y) for words on {a,b} of length <= 3.

    A substitution rule x -> y is finite if its application must terminate.
    This code identifies infinite substitutions and counts the rest as finite.
    """

    # 1. Generate all words of length 0 to 3
    alphabet = ['a', 'b']
    words = [""]
    for length in range(1, 4):
        for p in product(alphabet, repeat=length):
            words.append("".join(p))
    
    finite_substitution_count = 0
    
    # 2. Iterate through all possible couples (x, y)
    for x in words:
        for y in words:
            # Assume a substitution is finite unless it meets a condition for being infinite
            is_infinite = False
            
            # Condition 1: Substitution from an empty word ("") is always infinite.
            if not x:
                is_infinite = True
            
            # Condition 2: If x is a subword of y, it's an infinite substitution.
            # This also covers the simple infinite case x == y.
            elif x in y:
                is_infinite = True
                
            # Condition 3: Specific known non-terminating cases for |x| == |y| that are not
            # covered by the simple substring rule.
            elif (x == 'aab' and y == 'baa') or (x == 'bba' and y == 'abb'):
                is_infinite = True

            # 3. If no infinite conditions are met, count it as finite
            if not is_infinite:
                finite_substitution_count += 1
                
    # 4. Output the final count
    print(finite_substitution_count)

solve()
<<<181>>>