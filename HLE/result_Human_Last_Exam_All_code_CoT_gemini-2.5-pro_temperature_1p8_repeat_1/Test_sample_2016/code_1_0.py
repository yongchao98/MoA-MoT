def solve():
    """
    Simulates a plan to find the goal state.
    """
    # Define objects
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state
    state = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})"
    }

    # The plan to execute
    plan = [
        {'name': 'action_3', 'params': [o1, o2]},
        {'name': 'action_2', 'params': [o1, o1]},
        {'name': 'action_1', 'params': [o2, o1]},
        {'name': 'action_2', 'params': [o2, o2]},
        {'name': 'action_1', 'params': [o2, o1]},
        {'name': 'action_3', 'params': [o2, o2]},
        {'name': 'action_2', 'params': [o1, o2]},
        {'name': 'action_1', 'params': [o2, o2]},
        {'name': 'action_3', 'params': [o2, o1]},
        {'name': 'action_1', 'params': [o1, o2]},
        {'name': 'action_3', 'params': [o1, o1]},
    ]

    for step in plan:
        action_name = step['name']
        params = step['params']
        
        preconditions = set()
        add_effects = set()
        del_effects = set()

        if action_name == 'action_0':
            p0 = params[0]
            preconditions = {f"fluent_2({p0})"}
            del_effects = {f"fluent_2({p0})"}
        elif action_name == 'action_1':
            p0, p1 = params[0], params[1]
            # 'not fluent_2(p0)' is a negative precondition
            if f"fluent_2({p0})" in state:
                raise ValueError(f"Precondition 'not fluent_2({p0})' failed for {action_name}{params}")
            preconditions = {f"fluent_0({p0}, {p1})"}
            add_effects = {f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"}
            del_effects = {f"fluent_2({p1})", f"fluent_1({p1}, {p0})"}
        elif action_name == 'action_2':
            p0, p1 = params[0], params[1]
            preconditions = {f"fluent_2({p1})"}
            add_effects = {f"fluent_1({p0}, {p1})"}
            del_effects = {f"fluent_2({p1})"}
        elif action_name == 'action_3':
            p0, p1 = params[0], params[1]
            # 'not fluent_0(p1, p0)' is a negative precondition
            if f"fluent_0({p1}, {p0})" in state:
                raise ValueError(f"Precondition 'not fluent_0({p1}, {p0})' failed for {action_name}{params}")
            preconditions = {f"fluent_1({p1}, {p0})", f"fluent_2({p0})"}
            add_effects = {f"fluent_2({p0})", f"fluent_0({p1}, {p0})"}
            del_effects = {f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"}

        # Check positive preconditions
        if not preconditions.issubset(state):
            missing = preconditions - state
            raise ValueError(f"Precondition {missing} failed for {action_name}{params}")

        # Apply effects: deletes first, then adds
        state.difference_update(del_effects)
        state.update(add_effects)

    # Sort for consistent output and format as a string
    final_goal = "&".join(sorted(list(state)))
    print(final_goal)

solve()