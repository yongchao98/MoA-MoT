def solve():
    """
    Simulates the plan to determine the final goal state.
    """
    # Use a set to store the true fluents for efficient add/remove operations.
    current_state = {
        "fluent_1(object_1_type_0, object_2_type_0)",
        "fluent_1(object_2_type_0, object_1_type_0)",
        "fluent_2(object_1_type_0)",
    }

    # Define action effects
    def action_0(p0):
        # false effects
        current_state.discard(f"fluent_2({p0})")
        # true effects (none for this action)

    def action_1(p0, p1):
        # false effects
        current_state.discard(f"fluent_2({p1})")
        current_state.discard(f"fluent_1({p1}, {p0})")
        # true effects
        current_state.add(f"fluent_2({p0})")
        current_state.add(f"fluent_0({p0}, {p1})")
        current_state.add(f"fluent_1({p0}, {p1})")

    def action_2(p0, p1):
        # false effects
        current_state.discard(f"fluent_2({p1})")
        # true effects
        current_state.add(f"fluent_1({p0}, {p1})")

    def action_3(p0, p1):
        # false effects
        current_state.discard(f"fluent_0({p0}, {p1})")
        current_state.discard(f"fluent_1({p1}, {p0})")
        # true effects
        current_state.add(f"fluent_2({p0})")
        current_state.add(f"fluent_0({p1}, {p0})")

    # Define objects
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # The provided plan
    plan = [
        (action_3, (o1, o2)),
        (action_2, (o1, o1)),
        (action_1, (o2, o1)),
        (action_2, (o2, o2)),
        (action_1, (o2, o1)),
        (action_3, (o2, o2)),
        (action_2, (o1, o2)),
        (action_1, (o2, o2)),
        (action_3, (o2, o1)),
        (action_1, (o1, o2)),
        (action_3, (o1, o1)),
    ]

    # Execute the plan
    for action_func, params in plan:
        action_func(*params)

    # Sort the final fluents for a consistent order and format the output string
    goal_fluents = sorted(list(current_state))
    print("&".join(goal_fluents))

solve()