def solve():
    """
    Simulates a planning problem to determine the goal state.
    """

    def fluent_to_string(fluent_tuple):
        """Converts a fluent tuple to its string representation."""
        name, params = fluent_tuple
        return f"{name}({', '.join(params)})"

    def string_to_fluent(fluent_string):
        """Converts a fluent string back to a tuple representation."""
        name, params_str = fluent_string.split("(", 1)
        params = [p.strip() for p in params_str[:-1].split(",")]
        if params == ['']:
            params = []
        return (name, tuple(params))

    # Initial state represented as a set of true fluent strings
    initial_state_str = {
        "fluent_1(object_1_type_0, object_2_type_0)",
        "fluent_1(object_2_type_0, object_1_type_0)",
        "fluent_2(object_1_type_0)",
    }
    
    current_state = set(initial_state_str)

    # The sequence of actions in the plan
    plan = [
        ("action_3", ("object_1_type_0", "object_2_type_0")),
        ("action_2", ("object_1_type_0", "object_1_type_0")),
        ("action_1", ("object_2_type_0", "object_1_type_0")),
        ("action_2", ("object_2_type_0", "object_2_type_0")),
        ("action_1", ("object_2_type_0", "object_1_type_0")),
        ("action_3", ("object_2_type_0", "object_2_type_0")),
        ("action_2", ("object_1_type_0", "object_2_type_0")),
        ("action_1", ("object_2_type_0", "object_2_type_0")),
        ("action_3", ("object_2_type_0", "object_1_type_0")),
        ("action_1", ("object_1_type_0", "object_2_type_0")),
        ("action_3", ("object_1_type_0", "object_1_type_0")),
    ]

    # Process each action in the plan
    for action_name, params in plan:
        add_effects = set()
        del_effects = set()

        if action_name == "action_0":
            p0 = params[0]
            del_effects.add(f"fluent_2({p0})")
        
        elif action_name == "action_1":
            p0, p1 = params
            add_effects.add(f"fluent_2({p0})")
            add_effects.add(f"fluent_0({p0}, {p1})")
            add_effects.add(f"fluent_1({p0}, {p1})")
            del_effects.add(f"fluent_2({p1})")
            del_effects.add(f"fluent_1({p1}, {p0})")

        elif action_name == "action_2":
            p0, p1 = params
            add_effects.add(f"fluent_1({p0}, {p1})")
            del_effects.add(f"fluent_2({p1})")
            
        elif action_name == "action_3":
            p0, p1 = params
            add_effects.add(f"fluent_2({p0})")
            add_effects.add(f"fluent_0({p1}, {p0})")
            del_effects.add(f"fluent_0({p0}, {p1})")
            del_effects.add(f"fluent_1({p1}, {p0})")

        # Apply effects: first delete, then add
        current_state.difference_update(del_effects)
        current_state.update(add_effects)

    # Sort the final fluents for a consistent order
    sorted_goal_fluents = sorted(list(current_state))
    
    # Format the final goal string
    goal_string = "&".join(sorted_goal_fluents)
    print(goal_string)

solve()