def check_guess(guess, numbers, letters):
    # Check if a guess matches all the conditions
    def check_condition1(g):
        # 52FR: one number correct wrong pos, one too small, both letters wrong
        n_correct = sum(1 for i, n in enumerate(g[:2]) if n in numbers and n != numbers[i])
        n_small = sum(1 for n in g[:2] if int(n) < min(int(x) for x in numbers))
        return n_correct == 1 and n_small == 1 and g[2:] != letters

    def check_condition2(g):
        # 02IT: both numbers too small, both letters wrong
        return all(int(n) < min(int(x) for x in numbers) for n in g[:2])

    def check_condition3(g):
        # 37AO: both numbers wrong, both letters too early
        return (g[:2] != numbers) and all(l < min(letters) for l in g[2:])

    def check_condition4(g):
        # 17FJ: both numbers wrong, both letters too early
        return (g[:2] != numbers) and all(l < min(letters) for l in g[2:])

    def check_condition5(g):
        # 79NP: both numbers too large, one letter correct wrong pos
        return (all(int(n) > max(int(x) for x in numbers) for n in g[:2]) and
                sum(1 for i, l in enumerate(g[2:]) if l in letters and l != letters[i]) == 1)

    def check_condition6(g):
        # 91CK: both numbers wrong, both letters too early
        return (g[:2] != numbers) and all(l < min(letters) for l in g[2:])

    def check_condition7(g):
        # 56SM: both numbers correct wrong positions
        return set(g[:2]) == set(numbers) and g[:2] != numbers

    def check_condition8(g):
        # 76BJ: one number correct wrong pos, one too large
        n_correct = sum(1 for i, n in enumerate(g[:2]) if n in numbers and n != numbers[i])
        n_large = sum(1 for n in g[:2] if int(n) > max(int(x) for x in numbers))
        return n_correct == 1 and n_large == 1

    def check_condition9(g):
        # 31RG: both numbers too small, both letters wrong
        return all(int(n) < min(int(x) for x in numbers) for n in g[:2])

    def check_condition10(g):
        # 94AV: both numbers wrong, one letter correct right pos
        return (g[:2] != numbers and
                sum(1 for i, l in enumerate(g[2:]) if l == letters[i]) == 1)

    # Create the guess string
    guess_str = ''.join(map(str, guess))
    
    # Check all conditions
    return (check_condition1(guess_str) and check_condition2(guess_str) and
            check_condition3(guess_str) and check_condition4(guess_str) and
            check_condition5(guess_str) and check_condition6(guess_str) and
            check_condition7(guess_str) and check_condition8(guess_str) and
            check_condition9(guess_str) and check_condition10(guess_str))

# Generate all possible combinations and test them
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z
found = False

# From conditions 7, we know the numbers are 5 and 6
for n1 in ['5', '6']:
    for n2 in ['5', '6']:
        if n1 == n2:
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:
                    continue
                if check_guess([n1, n2, l1, l2], [n1, n2], [l1, l2]):
                    print([n1, n2, l1, l2])
                    found = True

if not found:
    print("No solution found")