def solve():
    """
    Simulates the plan from the initial state to find the goal state.
    """
    # Initial state as a set of true fluents
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)',
    }

    # The plan to execute
    plan = [
        "action_3(object_1_type_0, object_2_type_0)",
        "action_2(object_1_type_0, object_1_type_0)",
        "action_1(object_2_type_0, object_1_type_0)",
        "action_2(object_2_type_0, object_2_type_0)",
        "action_1(object_2_type_0, object_1_type_0)",
        "action_3(object_2_type_0, object_2_type_0)",
        "action_2(object_1_type_0, object_2_type_0)",
        "action_1(object_2_type_0, object_2_type_0)",
        "action_3(object_2_type_0, object_1_type_0)",
        "action_1(object_1_type_0, object_2_type_0)",
        "action_3(object_1_type_0, object_1_type_0)",
    ]

    for action_str in plan:
        # Parse action name and parameters
        action_name = action_str.split('(')[0]
        params = action_str.split('(')[1][:-1].split(', ')

        add_effects = []
        del_effects = []

        if action_name == "action_0":
            p0 = params[0]
            del_effects.append(f'fluent_2({p0})')
        
        elif action_name == "action_1":
            p0, p1 = params[0], params[1]
            add_effects.extend([
                f'fluent_2({p0})',
                f'fluent_0({p0}, {p1})',
                f'fluent_1({p0}, {p1})',
            ])
            del_effects.extend([
                f'fluent_2({p1})',
                f'fluent_1({p1}, {p0})',
            ])
            
        elif action_name == "action_2":
            p0, p1 = params[0], params[1]
            add_effects.append(f'fluent_1({p0}, {p1})')
            del_effects.append(f'fluent_2({p1})')

        elif action_name == "action_3":
            p0, p1 = params[0], params[1]
            add_effects.extend([
                f'fluent_2({p0})',
                f'fluent_0({p1}, {p0})',
            ])
            del_effects.extend([
                f'fluent_0({p0}, {p1})',
                f'fluent_1({p1}, {p0})',
            ])
        
        # Apply effects: process deletions first, then additions.
        for fluent in del_effects:
            state.discard(fluent)
        for fluent in add_effects:
            state.add(fluent)
            
    # Format the final state as the goal
    goal_fluents = sorted(list(state))
    print("&".join(goal_fluents))

solve()
<<<fluent_0(object_1_type_0, object_1_type_0)&fluent_0(object_1_type_0, object_2_type_0)&fluent_0(object_2_type_0, object_2_type_0)&fluent_1(object_1_type_0, object_2_type_0)&fluent_1(object_2_type_0, object_2_type_0)&fluent_2(object_1_type_0)>>>