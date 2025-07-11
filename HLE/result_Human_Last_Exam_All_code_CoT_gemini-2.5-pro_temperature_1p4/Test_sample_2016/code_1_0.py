def solve():
    """
    Simulates the plan to find the goal state.
    """
    # Define the objects for easier reference
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state
    current_state = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})",
    }

    # The optimal plan
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

    # Process each action in the plan
    for action_name, p0, p1 in plan:
        add_effects = set()
        del_effects = set()
        preconditions_met = False

        # Define action schemas
        if action_name == "action_0":
            if f"fluent_2({p0})" in current_state:
                preconditions_met = True
                del_effects.add(f"fluent_2({p0})")
        
        elif action_name == "action_1":
            if (f"fluent_2({p0})" not in current_state and 
                f"fluent_0({p0}, {p1})" in current_state):
                preconditions_met = True
                add_effects.add(f"fluent_2({p0})")
                add_effects.add(f"fluent_0({p0}, {p1})")
                add_effects.add(f"fluent_1({p0}, {p1})")
                del_effects.add(f"fluent_2({p1})")
                del_effects.add(f"fluent_1({p1}, {p0})")

        elif action_name == "action_2":
            if f"fluent_2({p1})" in current_state:
                preconditions_met = True
                add_effects.add(f"fluent_1({p0}, {p1})")
                del_effects.add(f"fluent_2({p1})")

        elif action_name == "action_3":
            if (f"fluent_1({p1}, {p0})" in current_state and 
                f"fluent_2({p0})" in current_state and 
                f"fluent_0({p1}, {p0})" not in current_state):
                preconditions_met = True
                add_effects.add(f"fluent_2({p0})")
                add_effects.add(f"fluent_0({p1}, {p0})")
                del_effects.add(f"fluent_0({p0}, {p1})")
                del_effects.add(f"fluent_1({p1}, {p0})")
        
        if not preconditions_met:
            # This block should not be reached if the plan is valid.
            print(f"Error: Preconditions not met for {action_name}({p0}, {p1})")
            return

        # Apply the effects: new_state = (current_state - deletes) U adds
        current_state = (current_state - del_effects).union(add_effects)

    # Sort the final fluents for a consistent output and join them with '&'
    goal_fluents = sorted(list(current_state))
    goal_string = "&".join(goal_fluents)
    print(goal_string)

solve()