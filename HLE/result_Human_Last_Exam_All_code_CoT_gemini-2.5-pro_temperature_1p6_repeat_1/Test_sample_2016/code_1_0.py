import re

def solve():
    """
    Simulates a planning problem to find the goal state.
    """
    # Initial state represented as a set of true fluents.
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)',
    }

    # The optimal plan to be executed.
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

    # Process each action in the plan
    for action_str in plan:
        # Parse the action string to get the name and parameters
        match = re.match(r'(\w+)\((.*)\)', action_str)
        action_name = match.group(1)
        params = [p.strip() for p in match.group(2).split(',')]

        add_effects = set()
        del_effects = set()

        # Define effects for each action type
        if action_name == 'action_0':
            p0 = params[0]
            del_effects.add(f'fluent_2({p0})')
        
        elif action_name == 'action_1':
            p0, p1 = params[0], params[1]
            add_effects.add(f'fluent_2({p0})')
            add_effects.add(f'fluent_0({p0}, {p1})')
            add_effects.add(f'fluent_1({p0}, {p1})')
            del_effects.add(f'fluent_2({p1})')
            del_effects.add(f'fluent_1({p1}, {p0})')

        elif action_name == 'action_2':
            p0, p1 = params[0], params[1]
            add_effects.add(f'fluent_1({p0}, {p1})')
            del_effects.add(f'fluent_2({p1})')

        elif action_name == 'action_3':
            p0, p1 = params[0], params[1]
            add_effects.add(f'fluent_2({p0})')
            add_effects.add(f'fluent_0({p1}, {p0})')
            del_effects.add(f'fluent_0({p0}, {p1})')
            del_effects.add(f'fluent_1({p1}, {p0})')

        # Apply the effects: delete effects first, then add effects.
        state.difference_update(del_effects)
        state.update(add_effects)

    # Sort the final fluents alphabetically and format the output string.
    goal_fluents = sorted(list(state))
    goal_str = '&'.join(goal_fluents)
    print(goal_str)

solve()