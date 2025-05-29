def check_conditions(password):
    # Convert password to list of characters
    pwd = list(password)
    
    # Condition 1: 87VG
    # - one number (7 or 8) is correct but in wrong position
    # - one number is incorrect and too large (must be 8)
    # - both letters are incorrect
    c1 = ('7' in pwd and pwd.index('7') != 0) and \
         '8' not in pwd and \
         'V' not in pwd and 'G' not in pwd

    # Condition 2: 53KJ
    # - both numbers are incorrect
    # - both letters are incorrect
    c2 = '5' not in pwd and '3' not in pwd and 'K' not in pwd and 'J' not in pwd

    # Condition 3: 68JZ
    # - both numbers are incorrect
    # - Z is correct and in correct position (position 3)
    # - J is incorrect
    c3 = '6' not in pwd and '8' not in pwd and \
         pwd[3] == 'Z' and 'J' not in pwd

    # Condition 4: 46JT
    # - both numbers are incorrect
    # - both letters are incorrect
    c4 = '4' not in pwd and '6' not in pwd and 'J' not in pwd and 'T' not in pwd

    # Condition 5: 39HZ
    # - both numbers are incorrect
    # - both letters (H and Z) are correct and in correct positions
    c5 = '3' not in pwd and '9' not in pwd and \
         pwd[2] == 'H' and pwd[3] == 'Z'

    # Condition 6: 83UC
    # - both numbers are incorrect
    # - both letters are incorrect
    c6 = '8' not in pwd and '3' not in pwd and 'U' not in pwd and 'C' not in pwd

    # Condition 7: 42CB
    # - one number is correct and in correct position
    # - one number is incorrect
    # - both letters are incorrect
    c7 = pwd[0] == '1' and '2' not in pwd and \
         'C' not in pwd and 'B' not in pwd

    return all([c1, c2, c3, c4, c5, c6, c7])

# We know:
# 1. H and Z are in positions 3 and 4 (from condition 5)
# 2. First number must be 1 (from condition 7)
# 3. 7 must be present but not in first position (from condition 1)
# 4. Many numbers are eliminated: 2,3,4,5,6,8,9

valid_numbers = '0123456789'
solutions = []

# Since we know first number is 1 and second must be 7
password = '17HZ'
if check_conditions(password):
    solutions.append(list(password))

print(solutions)