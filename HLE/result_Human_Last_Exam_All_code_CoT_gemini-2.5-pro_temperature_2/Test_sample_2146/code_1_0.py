import numpy as np
import itertools

def parse_state_string(s):
    """Parses the LaTeX-like matrix string into a numpy array."""
    s = s.strip()
    s = s.replace(r'\begin{array}{' + 'c'*20 + '}', '').replace(r'\end{array}', '')
    s = s.replace(r'\[', '').replace(r'\]', '')
    lines = s.strip().split(r'\\')
    grid = []
    for line in lines:
        line = line.strip()
        if not line: continue
        row = [int(x.strip()) for x in line.split('&')]
        grid.append(row)
    return np.array(grid, dtype=np.uint8)

states_str = {
1: r"0&0&0&1&1&1&1&1&0&0&0\\0&0&1&0&0&1&0&0&1&0&0\\0&1&1&1&1&0&1&1&1&1&0\\1&0&1&0&0&1&0&0&1&0&1\\1&0&1&0&0&1&0&0&1&0&1\\1&1&0&1&1&0&1&1&0&1&1\\1&0&1&0&0&1&0&0&1&0&1\\1&0&1&0&0&1&0&0&1&0&1\\0&1&1&1&1&0&1&1&1&1&0\\0&0&1&0&0&1&0&0&1&0&0\\0&0&0&1&1&1&1&1&0&0&0",
2: r"1&0&0&0&1&0&0&0&1&0&0&0&1\\0&1&1&1&0&0&0&0&0&1&1&1&0\\0&1&1&1&0&0&0&0&0&1&1&1&0\\0&1&1&1&1&1&1&1&1&1&1&1&0\\1&0&0&1&1&1&1&1&1&1&0&0&1\\0&0&0&1&1&1&1&1&1&1&0&0&0\\0&0&0&1&1&1&1&1&1&1&0&0&0\\0&0&0&1&1&1&1&1&1&1&0&0&0\\1&0&0&1&1&1&1&1&1&1&0&0&1\\0&1&1&1&1&1&1&1&1&1&1&1&0\\0&1&1&1&0&0&0&0&0&1&1&1&0\\0&1&1&1&0&0&0&0&0&1&1&1&0\\1&0&0&0&1&0&0&0&1&0&0&0&1",
3: r"0&0&1&0&0\\0&1&0&1&0\\1&0&1&0&1\\0&1&0&1&0\\0&0&1&0&0",
4: r"0&0&0&1&0&1&0&0&0\\0&1&1&0&0&0&1&1&0\\0&1&1&0&0&0&1&1&0\\1&0&0&0&0&0&0&0&1\\0&0&0&0&1&0&0&0&0\\1&0&0&0&0&0&0&0&1\\0&1&1&0&0&0&1&1&0\\0&1&1&0&0&0&1&1&0\\0&0&0&1&0&1&0&0&0",
5: r"1&1&1&0&0&1&0&0&1&1&1\\1&0&1&0&1&1&1&0&1&0&1\\1&1&1&0&0&1&0&0&1&1&1\\0&0&0&1&1&1&1&1&0&0&0\\0&1&0&1&1&0&1&1&0&1&0\\1&1&1&1&0&0&0&1&1&1&1\\0&1&0&1&1&0&1&1&0&1&0\\0&0&0&1&1&1&1&1&0&0&0\\1&1&1&0&0&1&0&0&1&1&1\\1&0&1&0&1&1&1&0&1&0&1\\1&1&1&0&0&1&0&0&1&1&1",
6: r"1&0&0&0&0&0&0&0&1\\0&1&1&0&0&0&1&1&0\\0&1&1&1&0&1&1&1&0\\0&0&1&1&1&1&1&0&0\\0&0&0&1&1&1&0&0&0\\0&0&1&1&1&1&1&0&0\\0&1&1&1&0&1&1&1&0\\0&1&1&0&0&0&1&1&0\\1&0&0&0&0&0&0&0&1",
7: r"1&1&0&0&0&1&1\\1&0&0&0&0&0&1\\0&0&1&1&1&0&0\\0&0&1&1&1&0&0\\0&0&1&1&1&0&0\\1&0&0&0&0&0&1\\1&1&0&0&0&1&1",
8: r"0&0&0&0&1&0&0&0&0\\0&0&1&1&0&1&1&0&0\\0&1&1&0&0&0&1&1&0\\0&1&0&1&1&1&0&1&0\\1&0&0&1&0&1&0&0&1\\0&1&0&1&1&1&0&1&0\\0&1&1&0&0&0&1&1&0\\0&0&1&1&0&1&1&0&0\\0&0&0&0&1&0&0&0&0",
9: r"1&0&1&0&1&0&1\\0&1&0&0&0&1&0\\1&0&0&1&0&0&1\\0&0&1&0&1&0&0\\1&0&0&1&0&0&1\\0&1&0&0&0&1&0\\1&0&1&0&1&0&1"
}

states = {k: parse_state_string(v.replace('c'*13, 'c'*20)) for k,v in states_str.items()}

def get_rule_constraints(prev_grid, next_grid):
    must_be_in = set()
    must_not_be_in = set()
    h1, w1 = prev_grid.shape
    h2, w2 = next_grid.shape
    pad_h = (h2 - h1) // 2 + 1
    pad_w = (w2 - w1) // 2 + 1
    padded_prev = np.pad(prev_grid, ((pad_h, pad_h), (pad_w, pad_w)), 'constant')
    h_pad, w_pad = padded_prev.shape
    start_r = (h_pad - h2) // 2
    start_c = (w_pad - w2) // 2

    for r_next in range(h2):
        for c_next in range(w2):
            r_prev_center = start_r + r_next
            c_prev_center = start_c + c_next
            neighborhood = padded_prev[r_prev_center-1 : r_prev_center+2, c_prev_center-1 : c_prev_center+2]
            s = np.sum(neighborhood)
            if next_grid[r_next, c_next] == 1: must_be_in.add(s)
            else: must_not_be_in.add(s)
    if not must_be_in.isdisjoint(must_not_be_in): return None
    return must_be_in, must_not_be_in

def combine_rules(rule1, rule2):
    in1, out1 = rule1
    in2, out2 = rule2
    new_in, new_out = in1.union(in2), out1.union(out2)
    if not new_in.isdisjoint(new_out): return None
    return new_in, new_out

def solve():
    s3_options, s4_options, s5_options = [7, 9], [4, 6, 8], [1, 5]
    r1_t3 = 7 
    s3_others = [s for s in s3_options if s != r1_t3]
    for r1_t4 in [4, 6, 8]: # Iterate all possibilities just in case my logic was flawed
        s4_others = [s for s in s4_options if s != r1_t4]
        for p_s4_others in itertools.permutations(s4_others):
            r2_t4, r3_t4 = p_s4_others
            for p_s5 in itertools.permutations(s5_options):
                r2_t5, r3_t5 = p_s5
                r1_seq = (3, r1_t3, r1_t4)
                r2_seq = (s3_others[0], r2_t4, r2_t5)
                r3_seq = (r3_t4, r3_t5, 2)
                sequences = {'r1': r1_seq, 'r2': r2_seq, 'r3': r3_seq}
                rules = {}
                valid_perm = True
                for name, seq in sequences.items():
                    part1 = get_rule_constraints(states[seq[0]], states[seq[1]])
                    if part1 is None: valid_perm = False; break
                    part2 = get_rule_constraints(states[seq[1]], states[seq[2]])
                    if part2 is None: valid_perm = False; break
                    rule = combine_rules(part1, part2)
                    if rule is None: valid_perm = False; break
                    rules[name] = rule
                if not valid_perm: continue
                if len(set(rules.values())) == 3:
                    rule1_labels = "".join(map(str, rules['r1'][1]))
                    rule2_labels = "".join(map(str, rules['r2'][1]))
                    rule3_labels = "".join(map(str, rules['r3'][1]))
                    print("{" + f"{''.join(map(str, r1_seq))},{''.join(map(str, r2_seq))},{''.join(map(str, r3_seq))}" + "}")
                    return

solve()