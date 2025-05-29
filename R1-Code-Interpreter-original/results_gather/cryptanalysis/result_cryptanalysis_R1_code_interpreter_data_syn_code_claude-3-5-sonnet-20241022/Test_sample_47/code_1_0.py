def check_conditions(password):
    # Convert password to list of characters
    pwd = list(password)
    
    # Condition 1: 87VG
    c1_part1 = ('7' in pwd and pwd.index('7') != 1) or ('8' in pwd and pwd.index('8') != 0)  # one number correct but wrong position
    c1_part2 = any(x in pwd for x in ['9'])  # one number too large
    c1_part3 = 'V' not in pwd and 'G' not in pwd  # both letters incorrect
    c1 = c1_part1 and c1_part2 and c1_part3

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