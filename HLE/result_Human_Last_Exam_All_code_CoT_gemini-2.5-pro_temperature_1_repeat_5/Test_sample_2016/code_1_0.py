def solve():
    """
    Simulates a planning problem to find the goal state.
    """

    # --- State Representation ---
    # The state is represented as a set of strings, where each string is a true fluent.
    state = {
        "fluent_1(object_1_type_0, object_2_type_0)",
        "fluent_1(object_2_type_0, object_1_type_0)",
        "fluent_2(object_1_type_0)",
    }

    # --- Object Definitions ---
    p = {
        "obj1": "object_1_type_0",
        "obj2": "object_2_type_0"
    }

    # --- Action Definitions ---

    def action_0(current_state, p0):
        # Preconditions: fluent_2(p0)
        if f"fluent_2({p0})" in current_state:
            new_state = current_state.copy()
            # Effects (del): fluent_2(p0)
            new_state.discard(f"fluent_2({p0})")
            return new_state
        # The plan is assumed to be valid, so no need for else.
        return current_state

    def action_1(current_state, p0, p1):
        # Preconditions: (not fluent_2(p0)) and fluent_0(p0, p1)
        if f"fluent_2({p0})" not in current_state and f"fluent_0({p0}, {p1})" in current_state:
            new_state = current_state.copy()
            # Effects (add): fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
            new_state.add(f"fluent_2({p0})")
            new_state.add(f"fluent_0({p0}, {p1})")
            new_state.add(f"fluent_1({p0}, {p1})")
            # Effects (del): fluent_2(p1), fluent_1(p1, p0)
            new_state.discard(f"fluent_2({p1})")
            new_state.discard(f"fluent_1({p1}, {p0})")
            return new_state
        return current_state

    def action_2(current_state, p0, p1):
        # Preconditions: fluent_2(p1)
        if f"fluent_2({p1})" in current_state:
            new_state = current_state.copy()
            # Effects (add): fluent_1(p0, p1)
            new_state.add(f"fluent_1({p0}, {p1})")
            # Effects (del): fluent_2(p1)
            new_state.discard(f"fluent_2({p1})")
            return new_state
        return current_state

    def action_3(current_state, p0, p1):
        # Preconditions: fluent_1(p1, p0) and fluent_2(p0) and (not fluent_0(p1, p0))
        if (f"fluent_1({p1}, {p0})" in current_state and
                f"fluent_2({p0})" in current_state and
                f"fluent_0({p1}, {p0})" not in current_state):
            new_state = current_state.copy()
            # Effects (add): fluent_2(p0), fluent_0(p1, p0)
            new_state.add(f"fluent_2({p0})")
            new_state.add(f"fluent_0({p1}, {p0})")
            # Effects (del): fluent_0(p0, p1), fluent_1(p1, p0)
            new_state.discard(f"fluent_0({p0}, {p1})")
            new_state.discard(f"fluent_1({p1}, {p0})")
            return new_state
        return current_state

    # --- Plan Execution ---
    plan = [
        lambda s: action_3(s, p["obj1"], p["obj2"]),
        lambda s: action_2(s, p["obj1"], p["obj1"]),
        lambda s: action_1(s, p["obj2"], p["obj1"]),
        lambda s: action_2(s, p["obj2"], p["obj2"]),
        lambda s: action_1(s, p["obj2"], p["obj1"]),
        lambda s: action_3(s, p["obj2"], p["obj2"]),
        lambda s: action_2(s, p["obj1"], p["obj2"]),
        lambda s: action_1(s, p["obj2"], p["obj2"]),
        lambda s: action_3(s, p["obj2"], p["obj1"]),
        lambda s: action_1(s, p["obj1"], p["obj2"]),
        lambda s: action_3(s, p["obj1"], p["obj1"]),
    ]

    final_state = state
    for action in plan:
        final_state = action(final_state)

    # --- Format and Print Goal ---
    goal_fluents = sorted(list(final_state))
    print("&".join(goal_fluents))

solve()