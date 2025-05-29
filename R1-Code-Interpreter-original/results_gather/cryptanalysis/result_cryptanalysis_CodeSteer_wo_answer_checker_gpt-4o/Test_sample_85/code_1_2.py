from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# All possible combinations of two numbers and two letters
possible_combinations = [''.join(p) for p in permutations(numbers, 2) for q in permutations(letters, 2)]

# Feedback constraints
constraints = [
    ('86AF', lambda n, l: n[0] not in '86' and n[1] not in '86' and l[0] not in 'AF' and l[1] not in 'AF'),
    ('98LF', lambda n, l: n[0] not in '98' and n[1] not in '98' and l[0] not in 'LF' and l[1] not in 'LF'),
    ('20XK', lambda n, l: (n[0] == '2' and n[1] != '0') or (n[0] != '2' and n[1] == '0') and l[0] not in 'XK' and l[1] not in 'XK'),
    ('38TM', lambda n, l: n[0] not in '38' and n[1] not in '38' and l[0] not in 'TM' and l[1] not in 'TM'),
    ('81JK', lambda n, l: n[0] not in '81' and n[1] not in '81' and l[0] not in 'JK' and l[1] not in 'JK'),
    ('23VQ', lambda n, l: (n[0] == '2' and n[1] != '3') or (n[0] != '2' and n[1] == '3') and l[0] not in 'VQ' and l[1] not in 'VQ'),
    ('85IC', lambda n, l: n[0] not in '85' and n[1] not in '85' and l[0] not in 'IC' and l[1] not in 'IC'),
    ('03OE', lambda n, l: n[0] not in '03' and n[1] not in '03' and l[0] not in 'OE' and l[1] not in 'OE'),
    ('78WK', lambda n, l: (n[0] in '78' and n[1] not in '78') or (n[0] not in '78' and n[1] in '78') and ((l[0] == 'W' and l[1] != 'K') or (l[0] != 'W' and l[1] == 'K'))),
    ('68KC', lambda n, l: n[0] not in '68' and n[1] not in '68' and l[0] not in 'KC' and l[1] not in 'KC'),
    ('47OM', lambda n, l: (n[0] == '7' and n[1] != '4') or (n[0] != '7' and n[1] == '4') and l[0] not in 'OM' and l[1] not in 'OM'),
    ('58ZQ', lambda n, l: n[0] not in '58' and n[1] not in '58' and l[0] not in 'ZQ' and l[1] not in 'ZQ'),
    ('49DW', lambda n, l: n[0] not in '49' and n[1] not in '49' and ((l[0] == 'W' and l[1] != 'D') or (l[0] != 'W' and l[1] == 'D'))),
]

# Filter combinations based on constraints
for guess, constraint in constraints:
    possible_combinations = [combo for combo in possible_combinations if constraint(combo[:2], combo[2:])]

# Output the valid combination
if possible_combinations:
    password = list(possible_combinations[0])
    print(f"<<< {password} >>>")
else:
    print("No valid combination found. Please check the constraints.")