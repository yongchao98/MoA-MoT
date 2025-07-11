def get_goal_from_plan():
    """
    This function simulates a plan to determine the final goal state.
    """

    # Initial state of the world
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)'
    }

    # Define the effects of each action
    def action_0(p0, current_state):
        to_add = set()
        to_remove = {f'fluent_2({p0})'}
        return (current_state - to_remove) | to_add

    def action_1(p0, p1, current_state):
        to_add = {
            f'fluent_2({p0})',
            f'fluent_0({p0}, {p1})',
            f'fluent_1({p0}, {p1})'
        }
        to_remove = {
            f'fluent_2({p1})',
            f'fluent_1({p1}, {p0})'
        }
        return (current_state - to_remove) | to_add

    def action_2(p0, p1, current_state):
        to_add = {f'fluent_1({p0}, {p1})'}
        to_remove = {f'fluent_2({p1})'}
        return (current_state - to_remove) | to_add

    def action_3(p0, p1, current_state):
        to_add = {
            f'fluent_2({p0})',
            f'fluent_0({p1}, {p0})'
        }
        to_remove = {
            f'fluent_0({p0}, {p1})',
            f'fluent_1({p1}, {p0})'
        }
        return (current_state - to_remove) | to_add

    # Execute the plan step by step
    state = action_3('object_1_type_0', 'object_2_type_0', state)
    state = action_2('object_1_type_0', 'object_1_type_0', state)
    state = action_1('object_2_type_0', 'object_1_type_0', state)
    state = action_2('object_2_type_0', 'object_2_type_0', state)
    state = action_1('object_2_type_0', 'object_1_type_0', state)
    state = action_3('object_2_type_0', 'object_2_type_0', state)
    state = action_2('object_1_type_0', 'object_2_type_0', state)
    state = action_1('object_2_type_0', 'object_2_type_0', state)
    state = action_3('object_2_type_0', 'object_1_type_0', state)
    state = action_1('object_1_type_0', 'object_2_type_0', state)
    state = action_3('object_1_type_0', 'object_1_type_0', state)

    # Sort the final fluents for consistent output
    goal_fluents = sorted(list(state))
    
    # Print the goal in the specified format
    print("&".join(goal_fluents))

get_goal_from_plan()
<<<fluent_0(object_1_type_0, object_1_type_0)&fluent_0(object_1_type_0, object_2_type_0)&fluent_0(object_2_type_0, object_2_type_0)&fluent_1(object_1_type_0, object_2_type_0)&fluent_1(object_2_type_0, object_2_type_0)&fluent_2(object_1_type_0)>>>