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
        f"fluent_2({o1})",
    }

    # The plan to execute
    plan = [
        ("action_3", o1, o2),
        ("action_2", o1, o1),
        ("action_1", o2, o1),
        ("action_2", o2, o2),
        ("action_1", o2, o1),
        ("action_3", o2, o2),
        ("action_2", o1, o2),
        ("action_1", o2, o2),
        ("action_3", o2, o1),
        ("action_1", o1, o2),
        ("action_3", o1, o1),
    ]

    def apply_action(current_state, action_name, params):
        """
        Applies a single action to the state.
        """
        preconditions = []
        add_effects = set()
        del_effects = set()

        p0 = params[0]
        if len(params) > 1:
            p1 = params[1]

        # Action Definitions
        if action_name == "action_0":
            # pre: fluent_2(p0)
            # del: fluent_2(p0)
            preconditions.append(f"fluent_2({p0})" in current_state)
            del_effects.add(f"fluent_2({p0})")

        elif action_name == "action_1":
            # pre: (not fluent_2(p0)), fluent_0(p0, p1)
            # add: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
            # del: fluent_2(p1), fluent_1(p1, p0)
            preconditions.append(f"fluent_2({p0})" not in current_state)
            preconditions.append(f"fluent_0({p0}, {p1})" in current_state)
            add_effects.add(f"fluent_2({p0})")
            add_effects.add(f"fluent_0({p0}, {p1})")
            add_effects.add(f"fluent_1({p0}, {p1})")
            del_effects.add(f"fluent_2({p1})")
            del_effects.add(f"fluent_1({p1}, {p0})")

        elif action_name == "action_2":
            # pre: fluent_2(p1)
            # add: fluent_1(p0, p1)
            # del: fluent_2(p1)
            preconditions.append(f"fluent_2({p1})" in current_state)
            add_effects.add(f"fluent_1({p0}, {p1})")
            del_effects.add(f"fluent_2({p1})")

        elif action_name == "action_3":
            # pre: fluent_1(p1, p0), fluent_2(p0), (not fluent_0(p1, p0))
            # add: fluent_2(p0), fluent_0(p1, p0)
            # del: fluent_0(p0, p1), fluent_1(p1, p0)
            preconditions.append(f"fluent_1({p1}, {p0})" in current_state)
            preconditions.append(f"fluent_2({p0})" in current_state)
            preconditions.append(f"fluent_0({p1}, {p0})" not in current_state)
            add_effects.add(f"fluent_2({p0})")
            add_effects.add(f"fluent_0({p1}, {p0})")
            del_effects.add(f"fluent_0({p0}, {p1})")
            del_effects.add(f"fluent_1({p1}, {p0})")

        # Check if all preconditions are met
        if not all(preconditions):
            raise ValueError(f"Precondition failed for {action_name}{params} in state {current_state}")

        # Apply effects: (state - dels) U adds
        new_state = (current_state - del_effects).union(add_effects)
        return new_state

    # Execute the plan
    for action in plan:
        action_name = action[0]
        params = action[1:]
        state = apply_action(state, action_name, params)
    
    # Format and print the final state (goal)
    goal_fluents = sorted(list(state))
    print("&".join(goal_fluents))

solve()