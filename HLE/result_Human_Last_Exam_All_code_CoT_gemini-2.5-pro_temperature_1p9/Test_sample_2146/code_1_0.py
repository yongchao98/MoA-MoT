import numpy as np
from scipy.signal import convolve2d

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by identifying the rule sequences based on
    grid size and transition consistency.
    """

    # Step 1: Define all state matrices from the problem.
    states = {
        1: np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]]),
        2: np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]]),
        3: np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]]),
        4: np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]]),
        5: np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]]),
        6: np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]]),
        7: np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]]),
        8: np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]]),
        9: np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]])
    }

    # Step 2: Implement a helper function to check sequence consistency.
    def find_sequence_consistency(seq):
        all_must_be, all_must_not_be = set(), set()
        kernel = np.ones((3, 3), dtype=int)
        for i in range(len(seq) - 1):
            grid_before, grid_after = states[seq[i]], states[seq[i+1]]
            if grid_after.shape[0] != grid_before.shape[0] + 2: return False
            
            sums = convolve2d(grid_before, kernel, mode='full')
            if sums.shape != grid_after.shape: return False
            
            must_be = set(sums[grid_after == 1])
            must_not_be = set(sums[grid_after == 0])
            
            if not must_be.isdisjoint(must_not_be): return False
            if not all_must_be.isdisjoint(must_not_be) or not all_must_not_be.isdisjoint(must_be): return False
            
            all_must_be.update(must_be)
            all_must_not_be.update(must_not_be)
            if not all_must_be.isdisjoint(all_must_not_be): return False
        return True

    # Step 3: Use size-to-time mapping to constrain the search.
    dim_to_time = {5: 2, 7: 3, 9: 4, 11: 5, 13: 6}
    L_t = {t: [] for t in range(2, 7)}
    for key, val in states.items():
        L_t[dim_to_time[val.shape[0]]].append(key)

    # Step 4: Systematically search for the sequences.
    # Find Rule 1: (t=2,3,4) -> (3, s3, s4)
    r1_seq = None
    for s3 in L_t[3]:
        for s4 in L_t[4]:
            seq = (3, s3, s4)
            if find_sequence_consistency(seq):
                r1_seq = seq
                break
        if r1_seq: break

    # Find Rule 3: (t=4,5,6) -> (s4, s5, 2) from remaining labels
    used_labels = set(r1_seq)
    r3_seq = None
    for s4 in [s for s in L_t[4] if s not in used_labels]:
        for s5 in [s for s in L_t[5] if s not in used_labels]:
            seq = (s4, s5, 2)
            if find_sequence_consistency(seq):
                r3_seq = seq
                break
        if r3_seq: break
        
    # Determine Rule 2 from the last remaining labels
    used_labels.update(r3_seq)
    rem_labels = set(range(1, 10)) - used_labels
    r2_t3 = [l for l in rem_labels if l in L_t[3]][0]
    r2_t4 = [l for l in rem_labels if l in L_t[4]][0]
    r2_t5 = [l for l in rem_labels if l in L_t[5]][0]
    r2_seq = (r2_t3, r2_t4, r2_t5)

    # Step 5: Format and print the final output
    print("Found matching sequences for the three rules:")
    print(f"Rule 1 (t=2,3,4) -> State labels: #{r1_seq[0]}, #{r1_seq[1]}, #{r1_seq[2]}")
    print(f"Rule 2 (t=3,4,5) -> State labels: #{r2_seq[0]}, #{r2_seq[1]}, #{r2_seq[2]}")
    print(f"Rule 3 (t=4,5,6) -> State labels: #{r3_seq[0]}, #{r3_seq[1]}, #{r3_seq[2]}")

    r1_str = "".join(map(str, r1_seq))
    r2_str = "".join(map(str, r2_seq))
    r3_str = "".join(map(str, r3_seq))

    final_answer = f"{{{r1_str},{r2_str},{r3_str}}}"
    print("\nFormatted Answer String:")
    print(final_answer)

solve_ca_puzzle()