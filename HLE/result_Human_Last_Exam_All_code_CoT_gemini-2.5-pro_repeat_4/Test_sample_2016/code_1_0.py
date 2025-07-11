def solve():
    """
    Simulates a planning problem to find the goal state.
    """

    # Define the objects
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state represented as a set of true fluents
    state = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})",
    }

    def apply_action(current_state, action_name, params):
        """
        Applies a single action to the state, checking preconditions and applying effects.
        """
        p = params
        new_state = set(current_state)

        if action_name == "action_0":
            # Preconditions for action_0(p0)
            if f"fluent_2({p[0]})" not in current_state:
                raise ValueError(f"Precondition failed for {action_name}{params}")
            # Effects
            new_state.discard(f"fluent_2({p[0]})")

        elif action_name == "action_1":
            # Preconditions for action_1(p0, p1)
            if (f"fluent_2({p[0]})" in current_state) or \
               (f"fluent_0({p[0]}, {p[1]})" not in current_state):
                raise ValueError(f"Precondition failed for {action_name}{params}")
            # Effects
            new_state.add(f"fluent_2({p[0]})")
            new_state.add(f"fluent_0({p[0]}, {p[1]})")
            new_state.add(f"fluent_1({p[0]}, {p[1]})")
            new_state.discard(f"fluent_2({p[1]})")
            new_state.discard(f"fluent_1({p[1]}, {p[0]})")

        elif action_name == "action_2":
            # Preconditions for action_2(p0, p1)
            if f"fluent_2({p[1]})" not in current_state:
                raise ValueError(f"Precondition failed for {action_name}{params}")
            # Effects
            new_state.add(f"fluent_1({p[0]}, {p[1]})")
            new_state.discard(f"fluent_2({p[1]})")

        elif action_name == "action_3":
            # Preconditions for action_3(p0, p1)
            if (f"fluent_1({p[1]}, {p[0]})" not in current_state) or \
               (f"fluent_2({p[0]})" not in current_state) or \
               (f"fluent_0({p[1]}, {p[0]})" in current_state):
                raise ValueError(f"Precondition failed for {action_name}{params}")
            # Effects
            new_state.add(f"fluent_2({p[0]})")
            new_state.add(f"fluent_0({p[1]}, {p[0]})")
            new_state.discard(f"fluent_0({p[0]}, {p[1]})")
            new_state.discard(f"fluent_1({p[1]}, {p[0]})")
            
        return new_state

    # The plan to execute
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

    # Execute the plan
    for action_name, params in plan:
        state = apply_action(state, action_name, params)

    # Sort the final fluents for a consistent output and format as a goal string
    goal_fluents = sorted(list(state))
    goal_string = "&".join(goal_fluents)
    print(goal_string)

solve()