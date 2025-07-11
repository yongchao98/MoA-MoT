def solve():
    """
    Simulates a planning problem to find the goal state.
    """
    # Initial state of the world
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)'
    }

    # Object names
    o1 = 'object_1_type_0'
    o2 = 'object_2_type_0'

    # The plan to execute
    plan = [
        {'name': 'action_3', 'params': (o1, o2)},
        {'name': 'action_2', 'params': (o1, o1)},
        {'name': 'action_1', 'params': (o2, o1)},
        {'name': 'action_2', 'params': (o2, o2)},
        {'name': 'action_1', 'params': (o2, o1)},
        {'name': 'action_3', 'params': (o2, o2)},
        {'name': 'action_2', 'params': (o1, o2)},
        {'name': 'action_1', 'params': (o2, o2)},
        {'name': 'action_3', 'params': (o2, o1)},
        {'name': 'action_1', 'params': (o1, o2)},
        {'name': 'action_3', 'params': (o1, o1)},
    ]

    for step in plan:
        action_name = step['name']
        params = step['params']
        
        preconditions_met = False
        adds = set()
        dels = set()

        if action_name == 'action_0':
            p0 = params[0]
            # Preconditions
            if f'fluent_2({p0})' in state:
                preconditions_met = True
                # Effects
                dels.add(f'fluent_2({p0})')

        elif action_name == 'action_1':
            p0, p1 = params
            # Preconditions
            if f'fluent_2({p0})' not in state and f'fluent_0({p0}, {p1})' in state:
                preconditions_met = True
                # Effects
                adds.add(f'fluent_2({p0})')
                adds.add(f'fluent_0({p0}, {p1})')
                adds.add(f'fluent_1({p0}, {p1})')
                dels.add(f'fluent_2({p1})')
                dels.add(f'fluent_1({p1}, {p0})')

        elif action_name == 'action_2':
            p0, p1 = params
            # Preconditions
            if f'fluent_2({p1})' in state:
                preconditions_met = True
                # Effects
                adds.add(f'fluent_1({p0}, {p1})')
                dels.add(f'fluent_2({p1})')
                
        elif action_name == 'action_3':
            p0, p1 = params
            # Preconditions
            if f'fluent_1({p1}, {p0})' in state and \
               f'fluent_2({p0})' in state and \
               f'fluent_0({p1}, {p0})' not in state:
                preconditions_met = True
                # Effects
                adds.add(f'fluent_2({p0})')
                adds.add(f'fluent_0({p1}, {p0})')
                dels.add(f'fluent_0({p0}, {p1})')
                dels.add(f'fluent_1({p1}, {p0})')

        if preconditions_met:
            state -= dels
            state |= adds
        else:
            # This case should not happen for a valid plan
            print(f"Error: Preconditions not met for {action_name}{params}")
            return

    # Sort the final state fluents alphabetically and print
    final_fluents = sorted(list(state))
    print("&".join(final_fluents))

solve()