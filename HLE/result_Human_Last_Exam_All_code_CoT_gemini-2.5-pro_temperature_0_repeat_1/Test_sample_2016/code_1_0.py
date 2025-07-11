def solve():
    """
    Simulates the plan to find the final goal state.
    """
    # Define the objects
    obj1 = "object_1_type_0"
    obj2 = "object_2_type_0"

    # Initial state represented as a set of true fluents
    state = {
        f"fluent_1({obj1}, {obj2})",
        f"fluent_1({obj2}, {obj1})",
        f"fluent_2({obj1})",
    }

    # The plan to execute
    plan = [
        ("action_3", (obj1, obj2)),
        ("action_2", (obj1, obj1)),
        ("action_1", (obj2, obj1)),
        ("action_2", (obj2, obj2)),
        ("action_1", (obj2, obj1)),
        ("action_3", (obj2, obj2)),
        ("action_2", (obj1, obj2)),
        ("action_1", (obj2, obj2)),
        ("action_3", (obj2, obj1)),
        ("action_1", (obj1, obj2)),
        ("action_3", (obj1, obj1)),
    ]

    # A helper function to apply actions
    def apply_action(current_state, action_name, params):
        p0 = params[0] if len(params) > 0 else None
        p1 = params[1] if len(params) > 1 else None
        
        preconditions_met = False
        add_effects = set()
        del_effects = set()

        if action_name == "action_0":
            # Preconditions: fluent_2(p0)
            if f"fluent_2({p0})" in current_state:
                preconditions_met = True
                # Effects: not fluent_2(p0)
                del_effects.add(f"fluent_2({p0})")

        elif action_name == "action_1":
            # Preconditions: (not fluent_2(p0)) and fluent_0(p0, p1)
            if f"fluent_2({p0})" not in current_state and f"fluent_0({p0}, {p1})" in current_state:
                preconditions_met = True
                # Effects: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1), not fluent_2(p1), not fluent_1(p1, p0)
                add_effects.add(f"fluent_2({p0})")
                add_effects.add(f"fluent_0({p0}, {p1})")
                add_effects.add(f"fluent_1({p0}, {p1})")
                del_effects.add(f"fluent_2({p1})")
                del_effects.add(f"fluent_1({p1}, {p0})")

        elif action_name == "action_2":
            # Preconditions: fluent_2(p1)
            if f"fluent_2({p1})" in current_state:
                preconditions_met = True
                # Effects: fluent_1(p0, p1), not fluent_2(p1)
                add_effects.add(f"fluent_1({p0}, {p1})")
                del_effects.add(f"fluent_2({p1})")

        elif action_name == "action_3":
            # Preconditions: fluent_1(p1, p0) and fluent_2(p0) and (not fluent_0(p1, p0))
            if (f"fluent_1({p1}, {p0})" in current_state and
                f"fluent_2({p0})" in current_state and
                f"fluent_0({p1}, {p0})" not in current_state):
                preconditions_met = True
                # Effects: fluent_2(p0), fluent_0(p1, p0), not fluent_0(p0, p1), not fluent_1(p1, p0)
                add_effects.add(f"fluent_2({p0})")
                add_effects.add(f"fluent_0({p1}, {p0})")
                del_effects.add(f"fluent_0({p0}, {p1})")
                del_effects.add(f"fluent_1({p1}, {p0})")

        if not preconditions_met:
            raise Exception(f"Precondition failed for {action_name}{params}")

        # Apply effects: remove del_effects then add add_effects
        new_state = (current_state - del_effects) | add_effects
        return new_state

    # Execute the plan
    for action_name, params in plan:
        state = apply_action(state, action_name, params)

    # Format the final state as the goal
    goal_fluents = sorted(list(state))
    goal_string = "&".join(goal_fluents)
    print(goal_string)

solve()