import numpy as np
from itertools import permutations

def solve():
    """
    Solves the cellular automaton puzzle by systematically searching for consistent
    state transition sequences for three distinct rules.
    """
    state_strings = {
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

    states = {}
    for label, text in state_strings.items():
        rows = text.strip().split('\n')
        parsed = [[int(c) for c in row.strip().split()] for row in rows]
        states[label] = np.array(parsed, dtype=np.int8)

    states_by_time = {
        2: [3],
        3: [7, 9],
        4: [4, 6, 8],
        5: [1, 5],
        6: [2]
    }

    def deduce_rule(prev_state, next_state):
        dim_prev = prev_state.shape[0]
        dim_next = next_state.shape[0]
        if dim_next != dim_prev + 2:
            return None

        on_set, off_set = set(), set()
        padded_prev = np.pad(prev_state, 1, 'constant', constant_values=0)
        grid_for_sum = np.pad(padded_prev, 1, 'constant', constant_values=0)

        for r in range(dim_next):
            for c in range(dim_next):
                s = np.sum(grid_for_sum[r:r + 3, c:c + 3])
                val = next_state[r, c]
                if val == 1:
                    on_set.add(s)
                else:
                    off_set.add(s)

        if not on_set.isdisjoint(off_set):
            return None
        return (on_set, off_set)

    def merge_rules(rule1, rule2):
        if rule1 is None or rule2 is None:
            return None
        on_set = rule1[0].union(rule2[0])
        off_set = rule1[1].union(rule2[1])
        if not on_set.isdisjoint(off_set):
            return None
        return (on_set, off_set)

    # Search for Rule 3: t=4,5,6
    r3_s6 = states_by_time[6][0]
    r3_s5, r3_s4, rule3 = None, None, None
    for s5_candidate in states_by_time[5]:
        rule_56 = deduce_rule(states[s5_candidate], states[r3_s6])
        if rule_56:
            for s4_candidate in states_by_time[4]:
                rule_45 = deduce_rule(states[s4_candidate], states[s5_candidate])
                if rule_45:
                    combined_rule = merge_rules(rule_45, rule_56)
                    if combined_rule:
                        r3_s5, r3_s4, rule3 = s5_candidate, s4_candidate, combined_rule
                        break
        if r3_s5:
            break
    
    rule3_labels = f"{r3_s4}{r3_s5}{r3_s6}"

    # Search for Rule 2: t=3,4,5
    rem_s5 = [s for s in states_by_time[5] if s != r3_s5]
    rem_s4 = [s for s in states_by_time[4] if s != r3_s4]
    
    r2_s5 = rem_s5[0]
    r2_s4, r2_s3, rule2 = None, None, None
    for s4_candidate in rem_s4:
        rule_45 = deduce_rule(states[s4_candidate], states[r2_s5])
        if rule_45:
            for s3_candidate in states_by_time[3]:
                rule_34 = deduce_rule(states[s3_candidate], states[s4_candidate])
                if rule_34:
                    combined_rule = merge_rules(rule_34, rule_45)
                    if combined_rule:
                        r2_s4, r2_s3, rule2 = s4_candidate, s3_candidate, combined_rule
                        break
        if r2_s4:
            break
            
    rule2_labels = f"{r2_s3}{r2_s4}{r2_s5}"

    # Deduce Rule 1: t=2,3,4
    r1_s2 = states_by_time[2][0]
    r1_s3 = [s for s in states_by_time[3] if s != r2_s3][0]
    r1_s4 = [s for s in rem_s4 if s != r2_s4][0]
    
    rule1_labels = f"{r1_s2}{r1_s3}{r1_s4}"
    
    # Final Answer construction
    final_answer = f"{{{rule1_labels},{rule2_labels},{rule3_labels}}}"
    print(final_answer)

solve()