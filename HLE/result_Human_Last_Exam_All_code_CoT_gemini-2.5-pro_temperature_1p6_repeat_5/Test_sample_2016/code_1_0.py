import collections

def get_goal_state():
    """
    Simulates the plan to determine the final goal state.
    """
    # Define object names for clarity
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state represented as a set of true fluents (facts)
    state = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})",
    }

    # The plan is a sequence of actions and their parameters
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
        # Create a copy of the state to modify
        next_state = state.copy()
        
        # Apply effects based on the action
        if action_name == "action_0":
            p0 = params[0]
            # Preconditions (for verification)
            assert f"fluent_2({p0})" in state
            # Delete effects
            next_state.discard(f"fluent_2({p0})")
            # Add effects: None

        elif action_name == "action_1":
            p0, p1 = params
            # Preconditions (for verification)
            assert f"fluent_2({p0})" not in state
            assert f"fluent_0({p0}, {p1})" in state
            # Delete effects
            next_state.discard(f"fluent_2({p1})")
            next_state.discard(f"fluent_1({p1}, {p0})")
            # Add effects
            next_state.add(f"fluent_2({p0})")
            next_state.add(f"fluent_0({p0}, {p1})")
            next_state.add(f"fluent_1({p0}, {p1})")
            
        elif action_name == "action_2":
            p0, p1 = params
            # Preconditions (for verification)
            assert f"fluent_2({p1})" in state
            # Delete effects
            next_state.discard(f"fluent_2({p1})")
            # Add effects
            next_state.add(f"fluent_1({p0}, {p1})")

        elif action_name == "action_3":
            p0, p1 = params
            # Preconditions (for verification)
            assert f"fluent_1({p1}, {p0})" in state
            assert f"fluent_2({p0})" in state
            assert f"fluent_0({p1}, {p0})" not in state
            # Delete effects
            next_state.discard(f"fluent_0({p0}, {p1})")
            next_state.discard(f"fluent_1({p1}, {p0})")
            # Add effects
            next_state.add(f"fluent_2({p0})")
            next_state.add(f"fluent_0({p1}, {p0})")

        # Update the state for the next step
        state = next_state

    # The goal is the final state. Sort fluents for a canonical representation.
    goal_fluents = sorted(list(state))
    
    # Format the output as a single string joined by '&'
    print("&".join(goal_fluents))

get_goal_state()