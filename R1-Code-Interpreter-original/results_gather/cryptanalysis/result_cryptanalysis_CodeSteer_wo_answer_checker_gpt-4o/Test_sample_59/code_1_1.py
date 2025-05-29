from itertools import permutations

# Define possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Define constraints based on feedback
constraints = [
    lambda x, y, z, w: x not in '20' and y not in '20' and z not in 'IX' and w not in 'IX',
    lambda x, y, z, w: (x == '0' or y == '0') and (x != '9' and y != '9') and z < 'Y' and w < 'T',
    lambda x, y, z, w: x not in '78' and y not in '78' and z > 'L' and w > 'K',
    lambda x, y, z, w: (x == '1' or y == '1') and (x != '5' and y != '5') and z > 'L' and w > 'D',
    lambda x, y, z, w: x not in '03' and y not in '03' and (z == 'C' or w == 'C') and (z > 'O' or w > 'O'),
    lambda x, y, z, w: (x == '1' or y == '1') and (x != '2' and y != '2') and z < 'X' and w < 'S',
    lambda x, y, z, w: x not in '48' and y not in '48' and z not in 'TP' and w not in 'TP',
    lambda x, y, z, w: x not in '25' and y not in '25' and z > 'L' and w > 'H',
    lambda x, y, z, w: x not in '72' and y not in '72' and (z == 'N' or w == 'N') and (z > 'O' or w > 'O'),
    lambda x, y, z, w: x not in '65' and y not in '65' and z < 'R' and w < 'X'
]

# Search for the correct combination
for x, y in permutations(numbers, 2):
    for z, w in permutations(letters, 2):
        if all(constraint(x, y, z, w) for constraint in constraints):
            print(f"<<< ['{x}', '{y}', '{z}', '{w}'] >>>")
            break