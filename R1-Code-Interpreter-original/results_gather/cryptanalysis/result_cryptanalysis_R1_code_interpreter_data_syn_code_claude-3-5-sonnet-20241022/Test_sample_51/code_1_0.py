def is_valid_solution(n1, n2, l1, l2):
    # Helper function to check if a letter comes before another in alphabet
    def is_before(a, b):
        return ord(a) < ord(b)

    # Condition 1: 52FR
    # - one number correct but wrong position (5)
    # - one number too small (2)
    # - both letters wrong
    c1 = ((n1 == '5' or n2 == '5') and int('2') < min(int(n1), int(n2)) and 
          l1 not in ['F', 'R'] and l2 not in ['F', 'R'])

    # Condition 2: 02IT
    # - both numbers too small
    c2 = int('0') < min(int(n1), int(n2)) and int('2') < min(int(n1), int(n2))

    # Condition 3: 37AO
    # - both numbers wrong
    # - both letters too early
    c3 = (n1 not in ['3', '7'] and n2 not in ['3', '7'] and 
          is_before('A', l1) and is_before('O', l2))

    # Condition 5: 79NP
    # - both numbers too large
    # - one letter correct but wrong position
    # - one letter too early
    c5 = (int('7') > max(int(n1), int(n2)) and int('9') > max(int(n1), int(n2)) and
          ((l1 == 'N' and l2 != 'N') or (l1 == 'P' and l2 != 'P') or
           (l2 == 'N' and l1 != 'N') or (l2 == 'P' and l1 != 'P')))

    # Condition 7: 56SM
    # - both numbers correct but wrong positions
    c7 = set([n1, n2]) == set(['5', '6'])

    # Condition 8: 76BJ
    # - one number correct but wrong position (6)
    # - one number too large (7)
    c8 = ((n1 == '6' or n2 == '6') and int('7') > max(int(n1), int(n2)))

    # Condition 10: 94AV
    # - both numbers wrong
    # - one letter correct and in correct position
    # - one letter too early
    c10 = (n1 not in ['9', '4'] and n2 not in ['9', '4'] and
           ((l1 == 'V' and l2 != 'A') or (l2 == 'V' and l1 != 'A')))

    return all([c1, c2, c3, c5, c7, c8, c10])

# Test all possible combinations with numbers 5 and 6
numbers = ['5', '6']
letters = [chr(i) for i in range(ord('A'), ord('Z')+1)]

for n1 in numbers:
    for n2 in numbers:
        if n1 == n2:
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:
                    continue
                if is_valid_solution(n1, n2, l1, l2):
                    print([n1, n2, l1, l2])