def solve_planning_goal():
    """
    Simulates a planning problem to determine the goal state.
    """
    # Initial state as a set of true fluents
    state = {
        "fluent_1(object_1_type_0, object_2_type_0)",
        "fluent_1(object_2_type_0, object_1_type_0)",
        "fluent_2(object_1_type_0)"
    }

    # The optimal plan provided
    plan = [
        ("action_3", "object_1_type_0", "object_2_type_0"),
        ("action_2", "object_1_type_0", "object_1_type_0"),
        ("action_1", "object_2_type_0", "object_1_type_0"),
        ("action_2", "object_2_type_0", "object_2_type_0"),
        ("action_1", "object_2_type_0", "object_1_type_0"),
        ("action_3", "object_2_type_0", "object_2_type_0"),
        ("action_2", "object_1_type_0", "object_2_type_0"),
        ("action_1", "object_2_type_0", "object_2_type_0"),
        ("action_3", "object_2_type_0", "object_1_type_0"),
        ("action_1", "object_1_type_0", "object_2_type_0"),
        ("action_3", "object_1_type_0", "object_1_type_0"),
    ]

    # Sequentially apply the effects of each action in the plan.
    # We don't need to check preconditions as the plan is given as valid.
    for action, p0, p1 in plan:
        if action == "action_0":
            # Effects-False: fluent_2(p0)
            state.discard(f"fluent_2({p0})")
        elif action == "action_1":
            # Effects-True: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
            # Effects-False: fluent_2(p1), fluent_1(p1, p0)
            state.add(f"fluent_2({p0})")
            state.add(f"fluent_0({p0}, {p1})")
            state.add(f"fluent_1({p0}, {p1})")
            state.discard(f"fluent_2({p1})")
            state.discard(f"fluent_1({p1}, {p0})")
        elif action == "action_2":
            # Effects-True: fluent_1(p0, p1)
            # Effects-False: fluent_2(p1)
            state.add(f"fluent_1({p0}, {p1})")
            state.discard(f"fluent_2({p1})")
        elif action == "action_3":
            # Effects-True: fluent_2(p0), fluent_0(p1, p0)
            # Effects-False: fluent_0(p0, p1), fluent_1(p1, p0)
            state.add(f"fluent_2({p0})")
            state.add(f"fluent_0({p1}, {p0})")
            state.discard(f"fluent_0({p0}, {p1})")
            state.discard(f"fluent_1({p1}, {p0})")

    # The goal is the final state. Sort fluents for a canonical representation.
    final_fluents = sorted(list(state))

    # Format the goal as a single string joined by '&'
    goal_string = "&".join(final_fluents)

    print(goal_string)

solve_planning_goal()