import numpy as np

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by identifying the evolution sequences for three distinct rules.
    """
    states = {
        1: np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]]),
        2: np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]]),
        3: np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]]),
        4: np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]]),
        5: np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]]),
        6: np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]]),
        7: np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]]),
        8: np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]]),
        9: np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]]),
    }

    def derive_rule(grid_A, grid_B):
        if grid_B.shape[0] != grid_A.shape[0] + 2: return None
        padded_A = np.pad(grid_A, 1, mode='constant', constant_values=0)
        rule = {}
        for r in range(grid_B.shape[0]):
            for c in range(grid_B.shape[1]):
                s = np.sum(padded_A[r:r+3, c:c+3])
                next_state = grid_B[r, c]
                if s in rule and rule[s] != next_state: return None
                rule[s] = next_state
        return rule

    def is_compatible(rule1, rule2):
        if not rule1 or not rule2: return False
        for s in rule1:
            if s in rule2 and rule1[s] != rule2[s]: return False
        return True

    t_groups = {
        2: [3], 3: [7, 9], 4: [4, 6, 8], 5: [1, 5], 6: [2]
    }
    results = {}

    # Find Rule 3 (t=4,5,6)
    found_rule3 = False
    for s4 in t_groups[4]:
        for s5 in t_groups[5]:
            s6 = t_groups[6][0]
            r1 = derive_rule(states[s4], states[s5])
            r2 = derive_rule(states[s5], states[s6])
            if is_compatible(r1, r2):
                results[3] = (s4, s5, s6)
                t_groups[4].remove(s4)
                t_groups[5].remove(s5)
                found_rule3 = True
                break
        if found_rule3: break

    # Find Rule 2 (t=3,4,5)
    found_rule2 = False
    for s3 in t_groups[3]:
        for s4 in t_groups[4]:
            s5 = t_groups[5][0]
            r1 = derive_rule(states[s3], states[s4])
            r2 = derive_rule(states[s4], states[s5])
            if is_compatible(r1, r2):
                results[2] = (s3, s4, s5)
                t_groups[3].remove(s3)
                t_groups[4].remove(s4)
                found_rule2 = True
                break
        if found_rule2: break
    
    # Rule 1 gets the remaining states
    results[1] = (t_groups[2][0], t_groups[3][0], t_groups[4][0])
    
    r1_labels = "".join(map(str, results[1]))
    r2_labels = "".join(map(str, results[2]))
    r3_labels = "".join(map(str, results[3]))

    print(f"{{{r1_labels},{r2_labels},{r3_labels}}")

solve_ca_puzzle()