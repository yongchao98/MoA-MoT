def solve():
    """
    Simulates the plan to find the goal state.
    """
    # Define the objects for clarity
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state as a set of true fluents
    state = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})",
    }

    # The plan of actions to execute
    plan = [
        ("action_3", (o1, o2)),
        ("action_2", (o1, o1)),
        ("action_1", (o2, o1)),
        ("action_2", (o2, o2)),
        ("action_1", (o2, o1)),
        ("action_3", (o2, o2)),
        ("action_2", (o1, o2)),
        ("action_1", (o2, o2)),
        ("action_3", (o2, o1)),
        ("action_1", (o1, o2)),
        ("action_3", (o1, o1)),
    ]

    for action_name, params in plan:
        add_effects = set()
        del_effects = set()

        if action_name == "action_0":
            p0, = params
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

        # Apply the effects: new_state = (current_state - del_effects) union add_effects
        state = (state - del_effects) | add_effects

    # Sort the final fluents for consistent output and join them
    goal_string = "&".join(sorted(list(state)))
    print(goal_string)

solve()