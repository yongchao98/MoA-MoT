from itertools import product

def test_combination(comb, conditions):
    n1, n2, l1, l2 = comb
    
    for guess, feedback in conditions:
        g_n1, g_n2 = int(guess[0]), int(guess[1])
        g_l1, g_l2 = guess[2], guess[3]
        
        # Check "86AF": both numbers incorrect, both letters incorrect
        if guess == "86AF":
            if n1 in (8,6) or n2 in (8,6) or l1 in ('A','F') or l2 in ('A','F'):
                return False
                
        # Check "98LF": numbers too large and incorrect, letters incorrect
        elif guess == "98LF":
            if n1 >= 9 or n2 >= 8 or l1 in ('L','F') or l2 in ('L','F'):
                return False
                
        # Check "20XK": one number correct and in position, one too small, letters incorrect
        elif guess == "20XK":
            correct_nums = sum(1 for x, y in zip([n1,n2], [2,0]) if x == y)
            if correct_nums != 1 or l1 in ('X','K') or l2 in ('X','K'):
                return False
                
        # Check "78WK": one number correct but wrong position, one letter correct in position
        elif guess == "78WK":
            if not ((n1 == 8 and n2 != 7) or (n2 == 7 and n1 != 8)):
                return False
            if not ((l1 == 'W' and l2 != 'K') or (l2 == 'K' and l1 != 'W')):
                return False
                
        # Check "47OM": one number correct in position
        elif guess == "47OM":
            correct_nums = sum(1 for x, y in zip([n1,n2], [4,7]) if x == y)
            if correct_nums != 1:
                return False
                
        # Check "49DW": numbers incorrect, one letter correct but wrong position
        elif guess == "49DW":
            if n1 in (4,9) or n2 in (4,9):
                return False
            if not ((l1 == 'W' and l2 != 'D') or (l2 == 'W' and l1 != 'D')):
                return False
    
    return True

# Generate possible combinations
numbers = range(10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
conditions = [
    ("86AF", "both incorrect"),
    ("98LF", "numbers too large"),
    ("20XK", "one number correct position"),
    ("78WK", "one number wrong position, one letter correct position"),
    ("47OM", "one number correct position"),
    ("49DW", "numbers incorrect, one letter wrong position")
]

solutions = []
for n1, n2 in product(range(8), range(8)):  # Numbers must be < 8 due to condition 2
    if n1 == n2:
        continue
    for l1, l2 in product(letters, letters):
        if l1 == l2:
            continue
        if test_combination((n1, n2, l1, l2), conditions):
            solutions.append([n1, n2, l1, l2])

print(f"Found {len(solutions)} solution(s):")
for sol in solutions:
    print(sol)