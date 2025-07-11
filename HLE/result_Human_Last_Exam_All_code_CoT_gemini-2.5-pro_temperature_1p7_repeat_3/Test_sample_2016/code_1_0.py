import collections

def solve():
    """
    Simulates a planning problem to determine the final goal state.
    """
    # Using defaultdict for convenience. A non-existent key will return an empty set.
    # The state is a dictionary mapping a fluent name to a set of true propositions.
    # e.g., state['fluent_0'] = {('obj1', 'obj2'), ...}
    state = collections.defaultdict(set)

    # Object constants
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"
    
    # --- Initial Conditions ---
    state['fluent_1'].add((o1, o2))
    state['fluent_1'].add((o2, o1))
    state['fluent_2'].add(o1)
    
    # --- Plan Definition ---
    plan = [
        ("action_3", (o1, o2)),
        ("action_2", (o1, o1)),
        ("action_1", (o2, o1)),
        ("action_2", (o2, o2)),
        ("action_1", (o2, o1)),
        ("action_3", (o2, o2)),
        ("action_2", (o1, o2)),
        ("action_1", (o2, o2)),
        ("action_3", (o2, o1)),
        ("action_1", (o1, o2)),
        ("action_3", (o1, o1)),
    ]
    
    # --- Action Execution Loop ---
    for action_name, params in plan:
        # Unpack parameters
        p0 = params[0]
        p1 = params[1] if len(params) > 1 else None

        if action_name == "action_0":
            # Preconditions: fluent_2(p0)
            # We assume the plan is valid, so no need to check preconditions.
            # Effects: not fluent_2(p0)
            state['fluent_2'].discard(p0)
        
        elif action_name == "action_1":
            # Preconditions: (not fluent_2(p0)) and fluent_0(p0, p1)
            # Effects: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
            #          not fluent_2(p1), not fluent_1(p1, p0)
            state['fluent_2'].add(p0)
            state['fluent_0'].add((p0, p1))
            state['fluent_1'].add((p0, p1))
            state['fluent_2'].discard(p1)
            state['fluent_1'].discard((p1, p0))

        elif action_name == "action_2":
            # Preconditions: fluent_2(p1)
            # Effects: fluent_1(p0, p1), not fluent_2(p1)
            state['fluent_1'].add((p0, p1))
            state['fluent_2'].discard(p1)

        elif action_name == "action_3":
            # Preconditions: fluent_1(p1, p0) and fluent_2(p0) and (not fluent_0(p1, p0))
            # Effects: fluent_2(p0), fluent_0(p1, p0)
            #          not fluent_0(p0, p1), not fluent_1(p1, p0)
            state['fluent_2'].add(p0)
            state['fluent_0'].add((p1, p0))
            state['fluent_0'].discard((p0, p1))
            state['fluent_1'].discard((p1, p0))
            
    # --- Format Final State ---
    final_facts = []
    # Sort for a canonical representation
    for p1, p2 in sorted(list(state['fluent_0'])):
        final_facts.append(f"fluent_0({p1}, {p2})")
    for p1, p2 in sorted(list(state['fluent_1'])):
        final_facts.append(f"fluent_1({p1}, {p2})")
    for p in sorted(list(state['fluent_2'])):
        final_facts.append(f"fluent_2({p})")

    # Print the final goal string
    print("&".join(final_facts))

solve()