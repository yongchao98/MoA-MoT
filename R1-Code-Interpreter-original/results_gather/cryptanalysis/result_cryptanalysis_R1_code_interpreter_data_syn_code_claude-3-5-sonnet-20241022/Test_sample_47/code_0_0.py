def check_conditions(password):
    # Convert password to list of characters
    pwd = list(password)
    
    # Condition 1: 87VG
    c1 = sum([
        (('8' in pwd or '7' in pwd) and pwd.index('8' if '8' in pwd else '7') != ('8' if '8' in pwd else '7' in pwd)),  # one number correct but wrong position
        any(n in pwd for n > '8'),  # one number too large
        'V' not in pwd and 'G' not in pwd  # both letters incorrect
    ]) == 3

    # Condition 2: 53KJ
    c2 = '5' not in pwd and '3' not in pwd and 'K' not in pwd and 'J' not in pwd

    # Condition 3: 68JZ
    c3 = '6' not in pwd and '8' not in pwd and \
         ((('J' in pwd and pwd.index('J') == 2) or ('Z' in pwd and pwd.index('Z') == 3)) and \
         not ('J' in pwd and 'Z' in pwd))

    # Condition 4: 46JT
    c4 = '4' not in pwd and '6' not in pwd and 'J' not in pwd and 'T' not in pwd

    # Condition 5: 39HZ
    c5 = '3' not in pwd and '9' not in pwd and \
         'H' in pwd and pwd.index('H') == 2 and \
         'Z' in pwd and pwd.index('Z') == 3

    # Condition 6: 83UC
    c6 = '8' not in pwd and '3' not in pwd and 'U' not in pwd and 'C' not in pwd

    # Condition 7: 42CB
    c7 = ((('4' in pwd and pwd.index('4') == 0) or ('2' in pwd and pwd.index('2') == 1)) and \
         not ('4' in pwd and '2' in pwd)) and \
         'C' not in pwd and 'B' not in pwd

    return all([c1, c2, c3, c4, c5, c6, c7])

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
solutions = []

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        password = n1 + n2 + l1 + l2
                        if check_conditions(password):
                            solutions.append([n1, n2, l1, l2])

print(solutions)