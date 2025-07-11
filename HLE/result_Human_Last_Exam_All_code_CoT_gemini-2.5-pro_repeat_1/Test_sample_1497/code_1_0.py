import sys

def solve():
    """
    This function solves the problem by identifying the parameter group affected by each variation.

    1. sul ponticello: Affects harmonic content -> ν -> group ii
    2. bridge mute: Affects bridge resonance -> f2, a2 -> group iv. Adding mass lowers frequency -> f2 goes down.
    3. helium room: Affects air resonance -> f1 -> group iii
    4. on E string: Affects open string fundamental -> F -> group i
    """
    
    # The groups corresponding to variations (1), (2), (3), and (4)
    group_1 = "ii"
    group_2 = "iv"
    group_3 = "iii"
    group_4 = "i"
    
    # The direction of change for the last parameter of the group for variation (2)
    direction_2 = "down"
    
    # Construct the final answer string
    final_answer = f"{group_1},{group_2},{group_3},{group_4},{direction_2}"
    
    print(final_answer)

solve()