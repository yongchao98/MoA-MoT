import numpy as np

def apply_permutation(p, perm_list):
    """Applies a permutation p to a list."""
    new_list = list(perm_list)
    # In a permutation tuple (a, b, c), a->b, b->c, c->a
    # We apply it to indices (0-based)
    for i in range(len(p) - 1):
        new_list[p[i]-1] = perm_list[p[i+1]-1]
    new_list[p[-1]-1] = perm_list[p[0]-1]
    return tuple(new_list)

def get_total_permutation(generators, initial_strands=(1, 2, 3, 4, 5)):
    """
    Calculates the total permutation of a braid word.
    Braid words are composed right-to-left (standard mathematical convention).
    """
    perm_map = {
        's1': (1, 2), 's1_inv': (1, 2),
        's2': (2, 3), 's2_inv': (2, 3),
        's3': (3, 4), 's3_inv': (3, 4),
        's4': (4, 5), 's4_inv': (4, 5),
    }

    # Powers of generators have identity permutation
    current_perm = list(initial_strands)
    
    # Process from right to left
    for gen_name in reversed(generators):
        # The permutation for sigma_i^k is the same as sigma_i if k is odd,
        # and identity if k is even. Here, exponents are 2, 2, 1, -1.
        if gen_name == 's1^2' or gen_name == 's2^2':
            # These have identity permutation, so they don't change strand positions
            perm = ()
        elif gen_name == 's3':
            perm = perm_map['s3']
        elif gen_name == 's4^-1':
            perm = perm_map['s4_inv']

        if perm:
            # Apply permutation by swapping elements based on positions
            # Let pi be the permutation (a, b). It swaps items at positions a-1 and b-1.
            a, b = perm
            current_perm[a-1], current_perm[b-1] = current_perm[b-1], current_perm[a-1]
            
    # The final permutation maps initial position i to the strand name at position i.
    # To find cycles, we find where each initial strand i ends up.
    final_pos = [0] * len(initial_strands)
    for i in range(len(initial_strands)):
        final_pos[current_perm[i]-1] = i + 1

    return final_pos

def find_cycles(p):
    """Finds cycles in a permutation list."""
    n = len(p)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if not visited[i]:
            cycle = []
            j = i
            while not visited[j]:
                visited[j] = True
                cycle.append(j + 1)
                j = p[j] - 1
            cycles.append(tuple(cycle))
    return cycles

# Braid word: sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1
braid_word = ['s1^2', 's2^2', 's3', 's4^-1']
# The permutation sends the strand at initial position i to a final position pi(i).
# Example: pi = [2, 3, 1] means 1->2, 2->3, 3->1. Cycle is (1,2,3).
# We calculate pi by tracking where strand j goes.
# Initial state: strands 1,2,3,4,5 are at positions 1,2,3,4,5
strands_at_pos = [1,2,3,4,5]

# Apply s4^-1: swaps strands at pos 4 and 5
strands_at_pos[3], strands_at_pos[4] = strands_at_pos[4], strands_at_pos[3]
# Now strands are (1,2,3,5,4) at positions (1,2,3,4,5)

# Apply s3: swaps strands at pos 3 and 4
strands_at_pos[2], strands_at_pos[3] = strands_at_pos[3], strands_at_pos[2]
# Now strands are (1,2,5,3,4) at positions (1,2,3,4,5)

# s2^2 and s1^2 have identity permutations, so no more changes.
final_strands_at_pos = strands_at_pos
# The permutation maps initial strand i to its final position.
# Strand 1 is at pos 1. Strand 2 is at pos 2. Strand 3 is at pos 4. Strand 4 is at pos 5. Strand 5 is at pos 3.
# So, pi(1)=1, pi(2)=2, pi(3)=4, pi(4)=5, pi(5)=3.
final_permutation = [0] * 5
for i in range(5):
    final_permutation[i] = final_strands_at_pos.index(i+1) + 1

cycles = find_cycles(final_permutation)

print("The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.")
print("The permutation induced by the braid determines the components of the link.")
print(f"The final permutation mapping initial position i to final position pi(i) is: {final_permutation}")
print(f"The cycle decomposition of the permutation is: {cycles}")
print("\nThis means the closure of the braid has 3 components:")
print(f"Component 1 is formed by strand {cycles[0][0]}.")
print(f"Component 2 is formed by strand {cycles[1][0]}.")
print(f"Component 3 is formed by strands {cycles[2]}.")

print("\nThe problem states two components are unknots. These are Components 1 and 2.")
print("Component 1 (strand 1) and Component 2 (strand 2) are simple loops linked with their neighbours, but have no self-crossings, so they are unknots.")

print("\nWe need to identify the knot type of Component 3, formed by strands 3, 4, and 5.")
print("The knot type is determined by the sub-braid acting on these three strands.")
print("The relevant parts of the braid word are sigma_3 and sigma_4^-1.")
print("Ignoring the linking with Component 2 (from sigma_2^2), we analyze the 3-braid formed by mapping strands {3,4,5} to {1,2,3}.")
print("Under this mapping, sigma_3 becomes sigma_1 and sigma_4^-1 becomes sigma_2^-1.")
print("The resulting 3-braid is sigma_1 * sigma_2^-1 (or sigma_2^-1 * sigma_1, which is conjugate and gives the same knot type).")

print("\nThe closure of the 3-braid sigma_1 * sigma_2^-1 is a well-known result in knot theory.")
print("It is the trefoil knot (3_1).")
print("\nTherefore, the other connected component is equivalent to a Trefoil knot.")
