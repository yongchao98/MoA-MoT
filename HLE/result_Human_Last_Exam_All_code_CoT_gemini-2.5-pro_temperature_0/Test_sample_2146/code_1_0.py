import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by finding a consistent assignment of states to rules.

    The method is based on a key insight: the standard "totalistic" rule as described
    leads to a logical contradiction. The puzzle is solvable under the assumption of an
    "outer totalistic" (or Life-like) rule, where the next state depends on the current
    cell's state and the sum of its 8 neighbors.

    This allows for a simplified analysis of the pattern's corners. Let c(t) be the
    state of the corner cell at time t. The corner at t+1 is born from a previously
    empty region, and its neighbor sum is simply c(t). This gives the relation:
    c(t+1) = 1 if c(t) is in the "Birth" set of the rule, and 0 otherwise.

    The script systematically checks all possible assignments of states to the timelines
    for each rule, and finds the unique assignment that is consistent with this
    corner-evolution logic for all three rules simultaneously.
    """

    # State properties: {label: (time, corner_value)}
    # Corner value is 1 if the corner cell of the matrix is 1, 0 otherwise.
    # Determined by inspecting the provided matrices.
    states = {
        1: (5, 0), 2: (6, 1), 3: (2, 0), 4: (4, 0),
        5: (5, 1), 6: (4, 1), 7: (3, 1), 8: (4, 0), 9: (3, 1)
    }

    # Group states by time
    states_by_time = {t: [] for t in range(2, 7)}
    for label, (time, corner) in states.items():
        states_by_time[time].append(label)

    # Get possible state sequences for each rule based on time
    r1_t2_options = states_by_time[2]
    r1_t3_options = states_by_time[3]
    r1_t4_options = states_by_time[4]

    r2_t3_options = states_by_time[3]
    r2_t4_options = states_by_time[4]
    r2_t5_options = states_by_time[5]

    r3_t4_options = states_by_time[4]
    r3_t5_options = states_by_time[5]
    r3_t6_options = states_by_time[6]

    # Generate all possible unique assignments of states to rules
    for r1_t2 in r1_t2_options:
        for r1_t3 in r1_t3_options:
            for r1_t4 in r1_t4_options:
                rule1_labels = {r1_t2, r1_t3, r1_t4}
                if len(rule1_labels) != 3: continue

                # Check Rule 1 consistency
                c2, c3, c4 = states[r1_t2][1], states[r1_t3][1], states[r1_t4][1]
                # c3 = (c2 in B1), c4 = (c3 in B1)
                b1_c2_is_1 = c3 == 1
                b1_c3_is_1 = c4 == 1
                if b1_c2_is_1 == (c2 == c3): continue # Contradiction: e.g. if c2=c3=1, then 1=(1 in B1) and 1=(1 in B1) is consistent. if c2=1,c3=0, then 0=(1 in B1).

                for r2_t3 in r2_t3_options:
                    for r2_t4 in r2_t4_options:
                        for r2_t5 in r2_t5_options:
                            rule2_labels = {r2_t3, r2_t4, r2_t5}
                            if len(rule2_labels) != 3: continue
                            if rule1_labels.intersection(rule2_labels): continue

                            # Check Rule 2 consistency
                            c3, c4, c5 = states[r2_t3][1], states[r2_t4][1], states[r2_t5][1]
                            # c4 = (c3 in B2), c5 = (c4 in B2)
                            b2_c3_is_1 = c4 == 1
                            b2_c4_is_1 = c5 == 1
                            if (c3 == c4 and b2_c3_is_1 != b2_c4_is_1): continue

                            # Remaining states for Rule 3
                            all_labels = set(states.keys())
                            rule3_labels_set = all_labels - rule1_labels - rule2_labels
                            
                            # Check if remaining states fit Rule 3's timeline
                            r3_t4, r3_t5, r3_t6 = -1, -1, -1
                            valid_r3 = True
                            for label in rule3_labels_set:
                                time = states[label][0]
                                if time == 4: r3_t4 = label
                                elif time == 5: r3_t5 = label
                                elif time == 6: r3_t6 = label
                                else: valid_r3 = False
                            if not valid_r3 or -1 in [r3_t4, r3_t5, r3_t6]: continue

                            # Check Rule 3 consistency
                            c4, c5, c6 = states[r3_t4][1], states[r3_t5][1], states[r3_t6][1]
                            # c5 = (c4 in B3), c6 = (c5 in B3)
                            b3_c4_is_1 = c5 == 1
                            b3_c5_is_1 = c6 == 1
                            if (c4 == c5 and b3_c4_is_1 != b3_c5_is_1): continue
                            
                            # Found a consistent assignment
                            r1_str = f"{r1_t2}{r1_t3}{r1_t4}"
                            r2_str = f"{r2_t3}{r2_t4}{r2_t5}"
                            r3_str = f"{r3_t4}{r3_t5}{r3_t6}"
                            
                            print(f"Found a consistent assignment:")
                            print(f"Rule 1 (t=2,3,4): States #{r1_t2}, #{r1_t3}, #{r1_t4}")
                            print(f"Rule 2 (t=3,4,5): States #{r2_t3}, #{r2_t4}, #{r2_t5}")
                            print(f"Rule 3 (t=4,5,6): States #{r3_t4}, #{r3_t5}, #{r3_t6}")
                            print("\nFinal Answer String:")
                            print(f"{{{r1_str},{r2_str},{r3_str}}}")
                            return

solve_ca_puzzle()