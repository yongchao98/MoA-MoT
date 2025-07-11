def solve():
    """
    Simulates a planning problem to find the goal state.
    """

    # Initial state represented as a set of true fluents
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)',
    }

    # The sequence of actions to be executed
    plan = [
        ('action_3', ('object_1_type_0', 'object_2_type_0')),
        ('action_2', ('object_1_type_0', 'object_1_type_0')),
        ('action_1', ('object_2_type_0', 'object_1_type_0')),
        ('action_2', ('object_2_type_0', 'object_2_type_0')),
        ('action_1', ('object_2_type_0', 'object_1_type_0')),
        ('action_3', ('object_2_type_0', 'object_2_type_0')),
        ('action_2', ('object_1_type_0', 'object_2_type_0')),
        ('action_1', ('object_2_type_0', 'object_2_type_0')),
        ('action_3', ('object_2_type_0', 'object_1_type_0')),
        ('action_1', ('object_1_type_0', 'object_2_type_0')),
        ('action_3', ('object_1_type_0', 'object_1_type_0')),
    ]

    # Definitions of actions with their preconditions and effects
    action_definitions = {
        'action_0': {
            'params': ['p0'],
            'pre': lambda s, p: f'fluent_2({p["p0"]})' in s,
            'add': lambda p: set(),
            'del': lambda p: {f'fluent_2({p["p0"]})'}
        },
        'action_1': {
            'params': ['p0', 'p1'],
            'pre': lambda s, p: (f'fluent_2({p["p0"]})' not in s and
                                 f'fluent_0({p["p0"]}, {p["p1"]})' in s),
            'add': lambda p: {f'fluent_2({p["p0"]})', f'fluent_0({p["p0"]}, {p["p1"]})', f'fluent_1({p["p0"]}, {p["p1"]})'},
            'del': lambda p: {f'fluent_2({p["p1"]})', f'fluent_1({p["p1"]}, {p["p0"]})'}
        },
        'action_2': {
            'params': ['p0', 'p1'],
            'pre': lambda s, p: f'fluent_2({p["p1"]})' in s,
            'add': lambda p: {f'fluent_1({p["p0"]}, {p["p1"]})'},
            'del': lambda p: {f'fluent_2({p["p1"]})'}
        },
        'action_3': {
            'params': ['p0', 'p1'],
            'pre': lambda s, p: (f'fluent_1({p["p1"]}, {p["p0"]})' in s and
                                 f'fluent_2({p["p0"]})' in s and
                                 f'fluent_0({p["p1"]}, {p["p0"]})' not in s),
            'add': lambda p: {f'fluent_2({p["p0"]})', f'fluent_0({p["p1"]}, {p["p0"]})'},
            'del': lambda p: {f'fluent_0({p["p0"]}, {p["p1"]})', f'fluent_1({p["p1"]}, {p["p0"]})'}
        }
    }

    # Execute the plan
    for action_name, action_params in plan:
        action = action_definitions[action_name]
        
        # Map parameters to their names
        param_map = {name: val for name, val in zip(action['params'], action_params)}

        # Check preconditions
        if not action['pre'](state, param_map):
            print(f"Error: Precondition failed for {action_name}{action_params}")
            return

        # Apply effects (delete then add)
        state -= action['del'](param_map)
        state |= action['add'](param_map)

    # Format the final goal state
    goal_fluents = sorted(list(state))
    goal_string = "&".join(goal_fluents)
    print(goal_string)

solve()