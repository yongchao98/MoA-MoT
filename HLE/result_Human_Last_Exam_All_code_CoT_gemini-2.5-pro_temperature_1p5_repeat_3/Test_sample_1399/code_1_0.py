import collections

def solve_cdcl_conflict():
    """
    Analyzes a specific CDCL conflict scenario to determine UIPs, the learned clause,
    and the backtrack level.
    """

    # --- 1. Model the System ---
    clauses = {
        1: [1, 6, 7],    # C1: x1 \/ x6 \/ x7
        2: [-2, -3, -4], # C2: not x2 \/ not x3 \/ not x4
        3: [5, -2],      # C3: x5 \/ not x2
        4: [4, -5, -6],  # C4: x4 \/ not x5 \/ not x6
        5: [6, -7],      # C5: x6 \/ not x7
    }
    decisions = [(-1, 1), (3, 2), (2, 3)] # (literal, level)

    # --- 2. Simulate BCP ---
    # assignments: var -> (value: bool, level: int, antecedent_clause_idx: int or None)
    assignments = {}
    # implication_order: tracks assignment order at the conflict level for resolution
    implication_order_l3 = []
    conflict_clause_idx = None
    conflict_level = 0
    q = collections.deque()

    for lit, level in decisions:
        current_level = level
        var = abs(lit)
        if var in assignments:
            continue
        
        assignments[var] = (lit > 0, level, None)
        q.append(lit)
        if level == 3:
            implication_order_l3.append(lit)

        # BCP loop
        while q:
            assigned_lit = q.popleft()
            
            for idx, clause in clauses.items():
                is_satisfied = any(
                    abs(l) in assignments and assignments[abs(l)][0] == (l > 0)
                    for l in clause
                )
                if is_satisfied:
                    continue

                unassigned_lits = [l for l in clause if abs(l) not in assignments]
                num_false_lits = sum(
                    1 for l in clause if abs(l) in assignments and assignments[abs(l)][0] != (l > 0)
                )

                if len(unassigned_lits) == 1:
                    unit_lit = unassigned_lits[0]
                    unit_var = abs(unit_lit)
                    if unit_var not in assignments:
                        assignments[unit_var] = (unit_lit > 0, current_level, idx)
                        q.append(unit_lit)
                        if current_level == 3:
                            implication_order_l3.append(unit_lit)
                elif num_false_lits == len(clause):
                    conflict_clause_idx = idx
                    conflict_level = current_level
                    break
            if conflict_clause_idx:
                break
        if conflict_clause_idx:
            break

    # --- 3. & 4. Conflict Analysis and 1UIP Clause Generation ---
    resolvent = set(clauses[conflict_clause_idx])

    while True:
        # Identify the literal in the resolvent that was assigned last at the conflict level
        last_assigned_var = None
        for implied_lit in reversed(implication_order_l3):
            if implied_lit in resolvent or -implied_lit in resolvent:
                last_assigned_var = abs(implied_lit)
                break
        
        # Resolve with the antecedent
        antecedent_idx = assignments[last_assigned_var][2]
        antecedent = set(clauses[antecedent_idx])
        
        # Perform resolution
        resolvent.remove(last_assigned_var)
        resolvent.remove(-last_assigned_var) # literals are guaranteed to exist for resolution
        antecedent.remove(last_assigned_var)
        antecedent.remove(-last_assigned_var)
        resolvent.update(antecedent)
        
        # Check for 1UIP condition
        num_lits_at_conflict_level = sum(
            1 for lit in resolvent if assignments[abs(lit)][1] == conflict_level
        )
        if num_lits_at_conflict_level == 1:
            learned_clause = resolvent
            break
            
    # --- 5. Identify Final Results ---

    # First UIP
    first_uip_lit_in_clause = next(
        lit for lit in learned_clause if assignments[abs(lit)][1] == conflict_level
    )
    first_uip_var = abs(first_uip_lit_in_clause)
    first_uip_assignment_val = assignments[first_uip_var][0]
    # The UIP assignment is the one that falsifies the literal in the clause
    # e.g., learned clause has `x6`, assignment was `x6=false`. UIP is `not x6`.
    first_uip = first_uip_var if not first_uip_assignment_val else -first_uip_var
    first_uip_str = f"{'not ' if first_uip < 0 else ''}x{abs(first_uip)}@{conflict_level}"

    # All UIPs
    decision_lit_at_conflict_level = decisions[conflict_level - 1][0]
    decision_uip_str = f"{'not ' if decision_lit_at_conflict_level < 0 else ''}x{abs(decision_lit_at_conflict_level)}@{conflict_level}"
    # The set of UIPs includes the first UIP and the decision literal.
    all_uips = sorted([first_uip_str.replace(" ", ""), decision_uip_str.replace(" ", "")])
    all_uips_str = ",".join(all_uips)

    # Learned Clause
    learned_clause_str = " \/ ".join(sorted(
        [f"{'not ' if l < 0 else ''}x{abs(l)}" for l in learned_clause],
        key=lambda x: int(x.strip('not x'))
    ))
    
    # Backtracking Level
    other_lits_levels = [
        assignments[abs(l)][1] for l in learned_clause if abs(l) != abs(first_uip)
    ]
    backtracking_level = max(other_lits_levels) if other_lits_levels else 0

    # --- 6. Format and Print Output ---
    final_answer = f"{all_uips_str},{first_uip_str.replace(' ', '')},{learned_clause_str.replace(' ', '')},{backtracking_level}"
    print(final_answer)

solve_cdcl_conflict()