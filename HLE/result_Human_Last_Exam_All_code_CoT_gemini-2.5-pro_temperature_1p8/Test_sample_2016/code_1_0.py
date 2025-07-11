def solve():
    """
    Simulates the provided plan to find the goal state.
    """
    # Define the objects for clarity
    o1 = 'object_1_type_0'
    o2 = 'object_2_type_0'

    # Initial state represented as a set of true fluents
    state = {
        f'fluent_1({o1}, {o2})',
        f'fluent_1({o2}, {o1})',
        f'fluent_2({o1})'
    }

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

    # Process each action in the plan
    for action in plan:
        adds = set()
        dels = set()
        p0, p1 = None, None
        
        if len(action['params']) > 0:
            p0 = action['params'][0]
        if len(action['params']) > 1:
            p1 = action['params'][1]

        if action['name'] == 'action_0':
            # del: fluent_2(p0)
            dels.add(f'fluent_2({p0})')
        
        elif action['name'] == 'action_1':
            # add: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
            adds.add(f'fluent_2({p0})')
            adds.add(f'fluent_0({p0}, {p1})')
            adds.add(f'fluent_1({p0}, {p1})')
            # del: fluent_2(p1), fluent_1(p1, p0)
            dels.add(f'fluent_2({p1})')
            dels.add(f'fluent_1({p1}, {p0})')

        elif action['name'] == 'action_2':
            # add: fluent_1(p0, p1)
            adds.add(f'fluent_1({p0}, {p1})')
            # del: fluent_2(p1)
            dels.add(f'fluent_2({p1})')
            
        elif action['name'] == 'action_3':
            # add: fluent_2(p0), fluent_0(p1, p0)
            adds.add(f'fluent_2({p0})')
            adds.add(f'fluent_0({p1}, {p0})')
            # del: fluent_0(p0, p1), fluent_1(p1, p0)
            dels.add(f'fluent_0({p0}, {p1})')
            dels.add(f'fluent_1({p1}, {p0})')
        
        # Apply the effects to the state
        state = (state - dels) | adds

    # Sort the final fluents for a consistent order and print
    final_goal = sorted(list(state))
    print("&".join(final_goal))

solve()