import collections

def set_to_str(clause_set):
    """Converts a set of literals like {1, -2} to 'x1 v not x2'."""
    parts = []
    # Sort for consistent ordering, by variable number
    for lit in sorted(list(clause_set), key=abs):
        var = abs(lit)
        if lit > 0:
            parts.append(f"x{var}")
        else:
            parts.append(f"not x{var}")
    return " v ".join(parts)

def lit_to_str(lit_assignment, level):
    """Converts a literal assignment like -6 at level 3 to 'not x6@3'."""
    var = abs(lit_assignment)
    level = assignment_info[lit_assignment][0]
    if lit_assignment > 0:
        return f"x{var}@{level}"
    else:
        return f"not x{var}@{level}"

def count_literals_at_level(clause_set, level, assignment_info):
    """Counts how many literals in a clause were assigned at a specific level."""
    count = 0
    for lit in clause_set:
        var = abs(lit)
        # Find the assignment for this variable (positive or negative)
        assignment = var if var in assignment_info else -var
        if assignment_info[assignment][0] == level:
            count += 1
    return count

# 1. Initial problem setup
clauses = {
    1: {1, 6, 7},
    2: {-2, -3, -4},
    3: {5, -2},
    4: {4, -5, -6},
    5: {6, -7}
}
conflict_level = 3

# 2. Simulate the propagation based on the given decisions.
# This trace shows the sequence of assignments.
# Format: (literal, decision_level, antecedent_clause_idx or 'D' for decision)
assignment_trace = collections.OrderedDict()
assignment_trace[-1] = (1, 'D')      # x1 = false @ 1
assignment_trace[3]  = (2, 'D')      # x3 = true @ 2
assignment_trace[2]  = (3, 'D')      # x2 = true @ 3
assignment_trace[5]  = (3, 3)      # C3(5 v -2) + (x2=T) => x5=T
assignment_trace[-4] = (3, 2)      # C2(-2 v -3 v -4) + (x2=T, x3=T) => x4=F
assignment_trace[-6] = (3, 4)      # C4(4 v -5 v -6) + (x4=F, x5=T) => x6=F
assignment_trace[7]  = (3, 1)      # C1(1 v 6 v 7) + (x1=F, x6=F) => x7=T

# The `assignment_trace` is ordered by assignment time. We need to look up info by literal.
assignment_info = {lit: (level, reason) for lit, (level, reason) in assignment_trace.items()}

# The conflict happens at clause C5 (x6 v -x7) because x6=F and x7=T.
conflict_clause = clauses[5].copy()

# 3. Perform conflict analysis using resolution.
print("--- Conflict Analysis ---")
print(f"Conflict found at Level {conflict_level} with Clause C5: {set_to_str(conflict_clause)}")
print("Assignments: x6 is false, x7 is true -> F v F\n")

# Start resolution with the conflict clause.
resolvent = conflict_clause
# Get a list of assignments made at the conflict level, in reverse chronological order.
level3_assignments = [lit for lit, (level, r) in assignment_trace.items() if level == 3]

# Loop until the clause is a 1UIP clause (only one literal from the current level).
while count_literals_at_level(resolvent, conflict_level, assignment_info) > 1:
    # Find the most recently assigned variable in the current resolvent to resolve away.
    latest_lit = None
    for lit in reversed(level3_assignments):
        if abs(lit) in {abs(r) for r in resolvent}:
            latest_lit = lit
            break

    var_to_resolve = abs(latest_lit)
    reason_clause_idx = assignment_info[latest_lit][1]
    reason_clause = clauses[reason_clause_idx]

    print(f"Resolving on variable x{var_to_resolve} (latest assignment was {lit_to_str(latest_lit, conflict_level)})")
    print(f"  Current clause: {set_to_str(resolvent)}")
    print(f"  Reason clause (C{reason_clause_idx}): {set_to_str(reason_clause)}")

    # Perform resolution: (resolvent U reason_clause) - {var, -var}
    resolvent = (resolvent.union(reason_clause)) - {var_to_resolve, -var_to_resolve}
    print(f"  New resolvent: {set_to_str(resolvent)}\n")

learned_clause = resolvent
print(f"Analysis complete. Learned clause: {set_to_str(learned_clause)}")


# 4. Extract answers from the analysis.
# The First UIP is the single literal remaining in the clause from the conflict level.
# Its assignment is what is falsifying the clause.
first_uip = None
for lit in learned_clause:
    assignment = lit if lit in assignment_info else -lit
    if assignment_info[assignment][0] == conflict_level:
        # The learned clause `x1 v x6` is false under assignments `x1=F` and `x6=F`.
        # The UIP is the literal at the conflict level that caused this, which is `not x6@3`.
        first_uip = -assignment
        break

# The only UIP at level 3 is not x6@3. In the implication graph, a path to the conflict can
# bypass x7, so x7 is not a UIP.
all_uips_str = lit_to_str(first_uip, conflict_level)
first_uip_str = lit_to_str(first_uip, conflict_level)
learned_clause_str = set_to_str(learned_clause)

# Backtracking level is the second-highest level in the learned clause.
levels_in_clause = [assignment_info[l if l in assignment_info else -l][0] for l in learned_clause]
levels_in_clause.sort(reverse=True)
backtrack_level = levels_in_clause[1] if len(levels_in_clause) > 1 else 0

final_answer = f"{all_uips_str}, {first_uip_str}, {learned_clause_str}, {backtrack_level}"
print("\n--- Final Answer ---")
print(final_answer)

# Final answer block
print(f"\n<<<{final_answer}>>>")