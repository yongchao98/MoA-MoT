def solve():
    """
    Simulates a plan to determine the final goal state.
    """
    # Initial state
    current_state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)',
    }

    # The provided optimal plan
    plan = [
        ('action_3', ['object_1_type_0', 'object_2_type_0']),
        ('action_2', ['object_1_type_0', 'object_1_type_0']),
        ('action_1', ['object_2_type_0', 'object_1_type_0']),
        ('action_2', ['object_2_type_0', 'object_2_type_0']),
        ('action_1', ['object_2_type_0', 'object_1_type_0']),
        ('action_3', ['object_2_type_0', 'object_2_type_0']),
        ('action_2', ['object_1_type_0', 'object_2_type_0']),
        ('action_1', ['object_2_type_0', 'object_2_type_0']),
        ('action_3', ['object_2_type_0', 'object_1_type_0']),
        ('action_1', ['object_1_type_0', 'object_2_type_0']),
        ('action_3', ['object_1_type_0', 'object_1_type_0']),
    ]

    for action_name, params in plan:
        add_effects = set()
        del_effects = set()

        if action_name == 'action_0':
            p0 = params[0]
            del_effects.add(f'fluent_2({p0})')

        elif action_name == 'action_1':
            p0, p1 = params
            add_effects.add(f'fluent_2({p0})')
            add_effects.add(f'fluent_0({p0}, {p1})')
            add_effects.add(f'fluent_1({p0}, {p1})')
            del_effects.add(f'fluent_2({p1})')
            del_effects.add(f'fluent_1({p1}, {p0})')

        elif action_name == 'action_2':
            p0, p1 = params
            add_effects.add(f'fluent_1({p0}, {p1})')
            del_effects.add(f'fluent_2({p1})')

        elif action_name == 'action_3':
            p0, p1 = params
            add_effects.add(f'fluent_2({p0})')
            add_effects.add(f'fluent_0({p1}, {p0})')
            del_effects.add(f'fluent_0({p0}, {p1})')
            del_effects.add(f'fluent_1({p1}, {p0})')
        
        # Apply effects: new_state = (current_state - del_effects) U add_effects
        current_state = (current_state - del_effects).union(add_effects)

    # Sort the final fluents for a consistent output format and join them
    goal_fluents = sorted(list(current_state))
    goal_string = "&".join(goal_fluents)
    print(goal_string)

solve()