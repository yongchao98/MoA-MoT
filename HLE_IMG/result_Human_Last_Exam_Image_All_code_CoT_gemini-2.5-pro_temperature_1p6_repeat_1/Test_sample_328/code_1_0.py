import collections

# There are 18 possible parameter configurations given the rules.
# This puzzle displays 16 of them.
# The unique parameter 'p' can be b, c, or d. Its value 'u' is in {-1, 0, 1}.
# The other two parameters have an equal value 'v' from {-1, 0, 1}, where v != u.

# Code format:
# B/b/z for unique b=1/-1/0
# C/c/Z for unique c=1/-1/0
# D/d/0 for unique d=1/-1/0

# This is a very challenging puzzle. The solution is found by careful deduction,
# identifying symmetric cases, inverse-color pairs, and characteristic behaviors
# like chaos, regular oscillation, and saturation.

# Let's define the parameters for each plot based on analysis.
# The plot number is the key (1-indexed), the value is the tuple (b,c,d)
params = {
    1: (0, -1, 0),   # Symmetric chaos -> pure c=-1
    2: (0, -1, -1),  # Blue-driven chaos
    3: (1, 0, 0),    # Pure red asymmetry
    4: (0, 1, 0),    # Symmetric regular oscillation -> pure c=1
    5: (0, 0, 1),    # Pure red saturation/drive
    6: (-1, 1, -1),  # Blue-driven regular oscillation
    7: (1, -1, 1),   # Red-driven chaos
    8: (1, 1, 0),    # Red-asymmetric regular oscillation
    9: (1, -1, -1),  # Blue-driven chaos with red asymmetry
    10: (-1, 0, 0),  # Pure blue asymmetry (color inverse of #3)
    11: (0, 1, 1),   # Red-driven regular oscillation
    12: (-1, 0, -1), # Strong blue saturation
    13: (-1, -1, 1), # Red-driven chaos with blue asymmetry
    14: (-1, -1, 0), # Blue-asymmetric chaos
    15: (0, 0, -1),  # Pure blue saturation/drive (color inverse of #5)
    16: (-1, 1, 1)   # Red-driven regular oscillation with blue asymmetry
}

# The two missing configs are (1,0,1) [Z] and (1,1,-1) [d].

def get_code(b, c, d):
    """Determines the code based on the unique parameter."""
    counts = collections.Counter([b, c, d])
    # Find the value that appears once (the unique one)
    unique_val = min(counts, key=counts.get)
    
    if b == unique_val and c == d:
        if b == 1: return 'B'
        if b == -1: return 'b'
        if b == 0: return 'z'
    elif c == unique_val and b == d:
        if c == 1: return 'C'
        if c == -1: return 'c'
        if c == 0: return 'Z'
    elif d == unique_val and b == c:
        if d == 1: return 'D'
        if d == -1: return 'd'
        if d == 0: return '0'
    return '?'

# Generate the 16-character answer string
answer_string = ""
for i in range(1, 17):
    b, c, d = params[i]
    answer_string += get_code(b, c, d)

print(answer_string)