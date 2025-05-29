from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Constraints based on feedback
constraints = [
    lambda x, y, z, w: x not in '87' and y not in '87' and z < 'V' and w < 'X',
    lambda x, y, z, w: x not in '26' and y not in '26' and z not in 'AP' and w not in 'AP',
    lambda x, y, z, w: (x == '9' and y not in '8') or (y == '8' and x not in '9') and z > 'Z' and w > 'T',
    lambda x, y, z, w: (x == '9' and y not in '8') or (y == '8' and x not in '9') and (z == 'H' or w == 'L'),
    lambda x, y, z, w: (x == '9' and y not in '3') or (y == '3' and x not in '9') and z not in 'BN' and w not in 'BN',
    lambda x, y, z, w: x not in '37' and y not in '37' and z > 'C' and w > 'G',
    lambda x, y, z, w: x not in '17' and y not in '17' and z > 'K' and w > 'W',
    lambda x, y, z, w: (x == '9' and y not in '2') or (y == '2' and x not in '9') and z not in 'XC' and w not in 'XC',
    lambda x, y, z, w: x not in '30' and y not in '30' and z < 'N' and w < 'X',
    lambda x, y, z, w: (x == '7' and y not in '4') or (y == '4' and x not in '7') and z not in 'DJ' and w not in 'DJ'
]

# Generate all possible combinations of numbers and letters
for x, y in permutations(numbers, 2):
    for z, w in permutations(letters, 2):
        if all(constraint(x, y, z, w) for constraint in constraints):
            print(f"<<< ['{x}', '{y}', '{z}', '{w}'] >>>")
            break