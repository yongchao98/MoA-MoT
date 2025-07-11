import numpy as np
from scipy.signal import convolve2d
from itertools import permutations

def parse_states():
    """Parses the state matrices from the problem description."""
    states_raw = {
        1: """
        0 0 0 1 1 1 1 1 0 0 0
        0 0 1 0 0 1 0 0 1 0 0
        0 1 1 1 1 0 1 1 1 1 0
        1 0 1 0 0 1 0 0 1 0 1
        1 0 1 0 0 1 0 0 1 0 1
        1 1 0 1 1 0 1 1 0 1 1
        1 0 1 0 0 1 0 0 1 0 1
        1 0 1 0 0 1 0 0 1 0 1
        0 1 1 1 1 0 1 1 1 1 0
        0 0 1 0 0 1 0 0 1 0 0
        0 0 0 1 1 1 1 1 0 0 0
        """,
        2: """
        1 0 0 0 1 0 0 0 1 0 0 0 1
        0 1 1 1 0 0 0 0 0 1 1 1 0
        0 1 1 1 0 0 0 0 0 1 1 1 0
        0 1 1 1 1 1 1 1 1 1 1 1 0
        1 0 0 1 1 1 1 1 1 1 0 0 1
        0 0 0 1 1 1 1 1 1 1 0 0 0
        0 0 0 1 1 1 1 1 1 1 0 0 0
        0 0 0 1 1 1 1 1 1 1 0 0 0
        1 0 0 1 1 1 1 1 1 1 0 0 1
        0 1 1 1 1 1 1 1 1 1 1 1 0
        0 1 1 1 0 0 0 0 0 1 1 1 0
        0 1 1 1 0 0 0 0 0 1 1 1 0
        1 0 0 0 1 0 0 0 1 0 0 0 1
        """,
        3: """
        0 0 1 0 0
        0 1 0 1 0
        1 0 1 0 1
        0 1 0 1 0
        0 0 1 0 0
        """,
        4: """
        0 0 0 1 0 1 0 0 0
        0 1 1 0 0 0 1 1 0
        0 1 1 0 0 0 1 1 0
        1 0 0 0 0 0 0 0 1
        0 0 0 0 1 0 0 0 0
        1 0 0 0 0 0 0 0 1
        0 1 1 0 0 0 1 1 0
        0 1 1 0 0 0 1 1 0
        0 0 0 1 0 1 0 0 0
        """,
        5: """
        1 1 1 0 0 1 0 0 1 1 1
        1 0 1 0 1 1 1 0 1 0 1
        1 1 1 0 0 1 0 0 1 1 1
        0 0 0 1 1 1 1 1 0 0 0
        0 1 0 1 1 0 1 1 0 1 0
        1 1 1 1 0 0 0 1 1 1 1
        0 1 0 1 1 0 1 1 0 1 0
        0 0 0 1 1 1 1 1 0 0 0
        1 1 1 0 0 1 0 0 1 1 1
        1 0 1 0 1 1 1 0 1 0 1
        1 1 1 0 0 1 0 0 1 1 1
        """,
        6: """
        1 0 0 0 0 0 0 0 1
        0 1 1 0 0 0 1 1 0
        0 1 1 1 0 1 1 1 0
        0 0 1 1 1 1 1 0 0
        0 0 0 1 1 1 0 0 0
        0 0 1 1 1 1 1 0 0
        0 1 1 1 0 1 1 1 0
        0 1 1 0 0 0 1 1 0
        1 0 0 0 0 0 0 0 1
        """,
        7: """
        1 1 0 0 0 1 1
        1 0 0 0 0 0 1
        0 0 1 1 1 0 0
        0 0 1 1 1 0 0
        0 0 1 1 1 0 0
        1 0 0 0 0 0 1
        1 1 0 0 0 1 1
        """,
        8: """
        0 0 0 0 1 0 0 0 0
        0 0 1 1 0 1 1 0 0
        0 1 1 0 0 0 1 1 0
        0 1 0 1 1 1 0 1 0
        1 0 0 1 0 1 0 0 1
        0 1 0 1 1 1 0 1 0
        0 1 1 0 0 0 1 1 0
        0 0 1 1 0 1 1 0 0
        0 0 0 0 1 0 0 0 0
        """,
        9: """
        1 0 1 0 1 0 1
        0 1 0 0 0 1 0
        1 0 0 1 0 0 1
        0 0 1 0 1 0 0
        1 0 0 1 0 0 1
        0 1 0 0 0 1 0
        1 0 1 0 1 0 1
        """
    }
    
    parsed = {}
    for i, s in states_raw.items():
        lines = s.strip().split('\n')
        parsed[i] = np.array([list(map(int, line.split())) for line in lines])
    return parsed

def trim(grid):
    """Trims rows and columns of zeros from a grid."""
    if not np.any(grid):
        return np.array([[0]])
    rows = np.any(grid, axis=1)
    cols = np.any(grid, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return grid[rmin:rmax+1, cmin:cmax+1]

def run_simulation(rule_tuple, max_t):
    """Simulates a CA rule starting from a single cell."""
    grid_size = 2 * max_t + 3
    grid = np.zeros((grid_size, grid_size), dtype=int)
    center = grid_size // 2
    grid[center, center] = 1
    
    kernel = np.ones((3, 3), dtype=int)
    
    patterns = {}
    for t in range(1, max_t + 1):
        sums = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
        grid = np.zeros_like(grid, dtype=int)
        for i in range(10):
            if rule_tuple[i] == 1:
                grid[sums == i] = 1
        
        trimmed_pattern = trim(grid)
        patterns[t] = trimmed_pattern
    return patterns

def main():
    target_states = parse_states()
    
    states_by_time = {}
    for label, matrix in target_states.items():
        t = (matrix.shape[0] - 1) // 2
        if t not in states_by_time:
            states_by_time[t] = []
        states_by_time[t].append((label, matrix))

    # Store rules that generate known patterns
    # rule_spec -> {t: label, ...}
    successful_rules = {}

    # R(0) must be 0 for finite patterns
    for i in range(1, 2**9):
        rule_code = i << 1
        rule_tuple = tuple((rule_code >> s) & 1 for s in range(10))
        
        # Simulate
        gen_patterns = run_simulation(rule_tuple, 6)
        
        # Match generated patterns to target states
        rule_matches = {}
        for t, pattern in gen_patterns.items():
            if t in states_by_time:
                for label, target_matrix in states_by_time[t]:
                    if pattern.shape == target_matrix.shape and np.array_equal(pattern, target_matrix):
                        rule_matches[t] = label
                        break
        
        if rule_matches:
            successful_rules[rule_tuple] = rule_matches

    # Find the combination of 3 distinct rules that partition the 9 labels
    
    # Possible sequences for Rule 1 (t=2,3,4)
    r1_candidates = []
    for rule, matches in successful_rules.items():
        if 2 in matches and 3 in matches and 4 in matches:
            r1_candidates.append({'rule': rule, 'seq': (matches[2], matches[3], matches[4])})

    # Possible sequences for Rule 2 (t=3,4,5)
    r2_candidates = []
    for rule, matches in successful_rules.items():
        if 3 in matches and 4 in matches and 5 in matches:
            r2_candidates.append({'rule': rule, 'seq': (matches[3], matches[4], matches[5])})

    # Possible sequences for Rule 3 (t=4,5,6)
    r3_candidates = []
    for rule, matches in successful_rules.items():
        if 4 in matches and 5 in matches and 6 in matches:
            r3_candidates.append({'rule': rule, 'seq': (matches[4], matches[5], matches[6])})

    # Brute-force the combination
    for r1_data in r1_candidates:
        for r2_data in r2_candidates:
            for r3_data in r3_candidates:
                # Rules must be distinct
                if r1_data['rule'] == r2_data['rule'] or \
                   r1_data['rule'] == r3_data['rule'] or \
                   r2_data['rule'] == r3_data['rule']:
                    continue

                seq1 = r1_data['seq']
                seq2 = r2_data['seq']
                seq3 = r3_data['seq']
                
                all_labels = list(seq1) + list(seq2) + list(seq3)
                
                # Labels must be a partition of {1..9}
                if len(set(all_labels)) == 9:
                    s1 = "".join(map(str, seq1))
                    s2 = "".join(map(str, seq2))
                    s3 = "".join(map(str, seq3))
                    
                    # Ensure first label of Rule 1 is #3 and last of Rule 3 is #2
                    if seq1[0] == 3 and seq3[2] == 2:
                        print(f"{{{s1},{s2},{s3}}}")
                        return

if __name__ == '__main__':
    main()