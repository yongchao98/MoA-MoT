import itertools

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite.
    A substitution is finite if any sequence of applications must terminate.
    Returns True if finite, False if infinite.
    """
    # Case 1: If x is the empty string, we can apply "" -> y infinitely.
    if not x:
        return False

    # Case 2: If the substitution makes the word shorter, it must terminate.
    if len(y) < len(x):
        return True

    # Case 3: If x -> x, it's an infinite loop on any word containing x.
    if x == y:
        return False
        
    # After the checks above, we are left with cases where len(y) >= len(x).
    # The substitution is infinite if it can "recreate" an occurrence of x.

    # Case 4: The substituted word y itself contains x.
    if x in y:
        return False

    # Case 5: A new x is created by an overlap. We test the minimal overlap
    # cases by checking the concatenated words 'x+y' and 'y+x'. If a new
    # copy of x is formed at the "seam", the number of occurrences of x
    # will be greater than the sum of occurrences in the component parts.
    
    # We check if (y+x) has more occurrences of x than y and x do individually.
    if (y + x).count(x) > y.count(x) + 1:
        return False
        
    # We check if (x+y) has more occurrences of x than x and y do individually.
    if (x + y).count(x) > 1 + y.count(x):
        return False

    # If none of the infinite conditions are met, the substitution is finite.
    return True

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3.
    """
    alphabet = ['a', 'b']
    words = ['']
    for length in range(1, 4):
        for item in itertools.product(alphabet, repeat=length):
            words.append("".join(item))

    # The number of such words is 1+2+4+8 = 15.
    # The number of pairs (x,y) is 15*15 = 225.
    # We will count how many of these 225 substitutions are finite.
    
    finite_count = 0
    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_count += 1
                
    print(f"{finite_count}")

solve()
<<<185>>>