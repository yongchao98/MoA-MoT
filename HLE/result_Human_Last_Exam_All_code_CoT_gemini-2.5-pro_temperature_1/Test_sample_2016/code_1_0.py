def solve():
    """
    Simulates a plan to determine the final goal state.
    """

    # Use a set to store the true fluents for efficient updates
    true_fluents = set()

    # Define objects
    obj1 = "object_1_type_0"
    obj2 = "object_2_type_0"

    # Set initial conditions
    true_fluents.add(f"fluent_1({obj1}, {obj2})")
    true_fluents.add(f"fluent_1({obj2}, {obj1})")
    true_fluents.add(f"fluent_2({obj1})")

    # --- Action definitions ---
    
    def action_0(p0):
        # Effects: del fluent_2(p0)
        true_fluents.discard(f"fluent_2({p0})")

    def action_1(p0, p1):
        # Effects (del): fluent_2(p1), fluent_1(p1, p0)
        true_fluents.discard(f"fluent_2({p1})")
        true_fluents.discard(f"fluent_1({p1}, {p0})")
        # Effects (add): fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
        true_fluents.add(f"fluent_2({p0})")
        true_fluents.add(f"fluent_0({p0}, {p1})")
        true_fluents.add(f"fluent_1({p0}, {p1})")

    def action_2(p0, p1):
        # Effects (del): fluent_2(p1)
        true_fluents.discard(f"fluent_2({p1})")
        # Effects (add): fluent_1(p0, p1)
        true_fluents.add(f"fluent_1({p0}, {p1})")

    def action_3(p0, p1):
        # Effects (del): fluent_0(p0, p1), fluent_1(p1, p0)
        true_fluents.discard(f"fluent_0({p0}, {p1})")
        true_fluents.discard(f"fluent_1({p1}, {p0})")
        # Effects (add): fluent_2(p0), fluent_0(p1, p0)
        true_fluents.add(f"fluent_2({p0})")
        true_fluents.add(f"fluent_0({p1}, {p0})")

    # --- Plan execution ---
    # The plan is assumed to be valid, so we don't check preconditions.

    action_3(obj1, obj2)
    action_2(obj1, obj1)
    action_1(obj2, obj1)
    action_2(obj2, obj2)
    action_1(obj2, obj1)
    action_3(obj2, obj2)
    action_2(obj1, obj2)
    action_1(obj2, obj2)
    action_3(obj2, obj1)
    action_1(obj1, obj2)
    action_3(obj1, obj1)

    # Sort the final fluents for a consistent output order
    final_goal_fluents = sorted(list(true_fluents))
    
    # Print the final goal in the specified format
    print("&".join(final_goal_fluents))

solve()
<<<fluent_0(object_1_type_0, object_1_type_0)&fluent_0(object_1_type_0, object_2_type_0)&fluent_0(object_2_type_0, object_2_type_0)&fluent_1(object_1_type_0, object_2_type_0)&fluent_1(object_2_type_0, object_2_type_0)&fluent_2(object_1_type_0)>>>