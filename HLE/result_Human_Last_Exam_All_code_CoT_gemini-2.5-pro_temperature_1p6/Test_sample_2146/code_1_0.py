import numpy as np
import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically checking all possible state transitions.
    """
    states = {
        1: np.array([[0,0,0,1,1,1,1,1,0,0,0], [0,0,1,0,0,1,0,0,1,0,0], [0,1,1,1,1,0,1,1,1,1,0], [1,0,1,0,0,1,0,0,1,0,1], [1,0,1,0,0,1,0,0,1,0,1], [1,1,0,1,1,0,1,1,0,1,1], [1,0,1,0,0,1,0,0,1,0,1], [1,0,1,0,0,1,0,0,1,0,1], [0,1,1,1,1,0,1,1,1,1,0], [0,0,1,0,0,1,0,0,1,0,0], [0,0,0,1,1,1,1,1,0,0,0]]),
        2: np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1], [0,1,1,1,0,0,0,0,0,1,1,1,0], [0,1,1,1,0,0,0,0,0,1,1,1,0], [0,1,1,1,1,1,1,1,1,1,1,1,0], [1,0,0,1,1,1,1,1,1,1,0,0,1], [0,0,0,1,1,1,1,1,1,1,0,0,0], [0,0,0,1,1,1,1,1,1,1,0,0,0], [0,0,0,1,1,1,1,1,1,1,0,0,0], [1,0,0,1,1,1,1,1,1,1,0,0,1], [0,1,1,1,1,1,1,1,1,1,1,1,0], [0,1,1,1,0,0,0,0,0,1,1,1,0], [0,1,1,1,0,0,0,0,0,1,1,1,0], [1,0,0,0,1,0,0,0,1,0,0,0,1]]),
        3: np.array([[0,0,1,0,0], [0,1,0,1,0], [1,0,1,0,1], [0,1,0,1,0], [0,0,1,0,0]]),
        4: np.array([[0,0,0,1,0,1,0,0,0], [0,1,1,0,0,0,1,1,0], [0,1,1,0,0,0,1,1,0], [1,0,0,0,0,0,0,0,1], [0,0,0,0,1,0,0,0,0], [1,0,0,0,0,0,0,0,1], [0,1,1,0,0,0,1,1,0], [0,1,1,0,0,0,1,1,0], [0,0,0,1,0,1,0,0,0]]),
        5: np.array([[1,1,1,0,0,1,0,0,1,1,1], [1,0,1,0,1,1,1,0,1,0,1], [1,1,1,0,0,1,0,0,1,1,1], [0,0,0,1,1,1,1,1,0,0,0], [0,1,0,1,1,0,1,1,0,1,0], [1,1,1,1,0,0,0,1,1,1,1], [0,1,0,1,1,0,1,1,0,1,0], [0,0,0,1,1,1,1,1,0,0,0], [1,1,1,0,0,1,0,0,1,1,1], [1,0,1,0,1,1,1,0,1,0,1], [1,1,1,0,0,1,0,0,1,1,1]]),
        6: np.array([[1,0,0,0,0,0,0,0,1], [0,1,1,0,0,0,1,1,0], [0,1,1,1,0,1,1,1,0], [0,0,1,1,1,1,1,0,0], [0,0,0,1,1,1,0,0,0], [0,0,1,1,1,1,1,0,0], [0,1,1,1,0,1,1,1,0], [0,1,1,0,0,0,1,1,0], [1,0,0,0,0,0,0,0,1]]),
        7: np.array([[1,1,0,0,0,1,1], [1,0,0,0,0,0,1], [0,0,1,1,1,0,0], [0,0,1,1,1,0,0], [0,0,1,1,1,0,0], [1,0,0,0,0,0,1], [1,1,0,0,0,1,1]]),
        8: np.array([[0,0,0,0,1,0,0,0,0], [0,0,1,1,0,1,1,0,0], [0,1,1,0,0,0,1,1,0], [0,1,0,1,1,1,0,1,0], [1,0,0,1,0,1,0,0,1], [0,1,0,1,1,1,0,1,0], [0,1,1,0,0,0,1,1,0], [0,0,1,1,0,1,1,0,0], [0,0,0,0,1,0,0,0,0]]),
        9: np.array([[1,0,1,0,1,0,1], [0,1,0,0,0,1,0], [1,0,0,1,0,0,1], [0,0,1,0,1,0,0], [1,0,0,1,0,0,1], [0,1,0,0,0,1,0], [1,0,1,0,1,0,1]])
    }

    by_time = {
        2: [3],
        3: [7, 9],
        4: [4, 6, 8],
        5: [1, 5],
        6: [2]
    }

    def check_transition(grid_t, grid_t1):
        if grid_t1.shape[0] != grid_t.shape[0] + 2:
            return None
        
        padded_t = np.pad(grid_t, ((1, 1), (1, 1)), 'constant', constant_values=0)
        rule = {}
        
        for r, c in np.ndindex(grid_t1.shape):
            sum_val = np.sum(padded_t[r:r+3, c:c+3])
            target_state = grid_t1[r, c]
            if sum_val in rule and rule[sum_val] != target_state:
                return None  # Inconsistent rule
            rule[sum_val] = target_state
        return rule

    def merge_rules(rule1, rule2):
        merged = rule1.copy()
        for s, v in rule2.items():
            if s in merged and merged[s] != v:
                return None  # Conflict
            merged[s] = v
        return merged

    # Find the state sequence for each rule
    
    # Rule 1: t=2,3,4
    s3_options = by_time[3]
    s4_options = by_time[4]
    
    rule1_seq = None
    for s3 in s3_options:
        for s4 in s4_options:
            rule_23 = check_transition(states[by_time[2][0]], states[s3])
            if rule_23 is None: continue
            
            rule_34 = check_transition(states[s3], states[s4])
            if rule_34 is None: continue
            
            if merge_rules(rule_23, rule_34) is not None:
                rule1_seq = (by_time[2][0], s3, s4)
                break
        if rule1_seq: break

    # Rule 2: t=3,4,5
    s3_options_r2 = [s for s in by_time[3] if s not in rule1_seq]
    s4_options_r2 = [s for s in by_time[4] if s not in rule1_seq]
    s5_options_r2 = by_time[5]
    
    rule2_seq = None
    for s4 in s4_options_r2:
        for s5 in s5_options_r2:
            rule_34 = check_transition(states[s3_options_r2[0]], states[s4])
            if rule_34 is None: continue

            rule_45 = check_transition(states[s4], states[s5])
            if rule_45 is None: continue

            if merge_rules(rule_34, rule_45) is not None:
                # To resolve ambiguity, we check Rule 3. A valid Rule 3 chain must exist.
                s4_options_r3 = [s for s in s4_options_r2 if s != s4]
                s5_options_r3 = [s for s in s5_options_r2 if s != s5]
                
                # Check if a valid Rule 3 can be formed with the remaining states
                r45_r3 = check_transition(states[s4_options_r3[0]], states[s5_options_r3[0]])
                if r45_r3 is None: continue
                r56_r3 = check_transition(states[s5_options_r3[0]], states[by_time[6][0]])
                if r56_r3 is None: continue
                if merge_rules(r45_r3, r56_r3) is not None:
                    rule2_seq = (s3_options_r2[0], s4, s5)
                    break
        if rule2_seq: break
        
    # Rule 3: t=4,5,6
    s4_options_r3 = [s for s in by_time[4] if s not in rule1_seq and s not in rule2_seq]
    s5_options_r3 = [s for s in by_time[5] if s not in rule2_seq]
    rule3_seq = (s4_options_r3[0], s5_options_r3[0], by_time[6][0])
    
    # Format the final answer
    r1_str = "".join(map(str, rule1_seq))
    r2_str = "".join(map(str, rule2_seq))
    r3_str = "".join(map(str, rule3_seq))

    print(f"{{{r1_str},{r2_str},{r3_str}}}")

solve_ca_puzzle()