def solve():
    """
    Simulates the given plan to find the goal state.
    """
    # Define object names for convenience
    obj1 = "object_1_type_0"
    obj2 = "object_2_type_0"

    # Initial state represented as a set of true fluents
    state = {
        f"fluent_1({obj1}, {obj2})",
        f"fluent_1({obj2}, {obj1})",
        f"fluent_2({obj1})"
    }

    # --- Helper function to apply actions ---
    def apply_action(preconditions_met, add_effects, del_effects):
        # The plan is optimal, so preconditions are assumed to be met.
        # This function applies the changes to the state.
        for fluent in del_effects:
            state.discard(fluent)
        for fluent in add_effects:
            state.add(fluent)

    # --- Plan Execution ---

    # 1. action_3(object_1_type_0, object_2_type_0)
    p0, p1 = obj1, obj2
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p1}, {p0})"},
        del_effects={f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"}
    )

    # 2. action_2(object_1_type_0, object_1_type_0)
    p0, p1 = obj1, obj1
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})"}
    )

    # 3. action_1(object_2_type_0, object_1_type_0)
    p0, p1 = obj2, obj1
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})", f"fluent_1({p1}, {p0})"}
    )

    # 4. action_2(object_2_type_0, object_2_type_0)
    p0, p1 = obj2, obj2
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})"}
    )

    # 5. action_1(object_2_type_0, object_1_type_0)
    p0, p1 = obj2, obj1
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})", f"fluent_1({p1}, {p0})"}
    )

    # 6. action_3(object_2_type_0, object_2_type_0)
    p0, p1 = obj2, obj2
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p1}, {p0})"},
        del_effects={f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"}
    )
    
    # 7. action_2(object_1_type_0, object_2_type_0)
    p0, p1 = obj1, obj2
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})"}
    )

    # 8. action_1(object_2_type_0, object_2_type_0)
    p0, p1 = obj2, obj2
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})", f"fluent_1({p1}, {p0})"}
    )

    # 9. action_3(object_2_type_0, object_1_type_0)
    p0, p1 = obj2, obj1
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p1}, {p0})"},
        del_effects={f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"}
    )

    # 10. action_1(object_1_type_0, object_2_type_0)
    p0, p1 = obj1, obj2
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"},
        del_effects={f"fluent_2({p1})", f"fluent_1({p1}, {p0})"}
    )

    # 11. action_3(object_1_type_0, object_1_type_0)
    p0, p1 = obj1, obj1
    apply_action(
        preconditions_met=True,
        add_effects={f"fluent_2({p0})", f"fluent_0({p1}, {p0})"},
        del_effects={f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"}
    )

    # Sort the final fluents for a canonical goal representation and print
    goal_fluents = sorted(list(state))
    print("&".join(goal_fluents))

solve()