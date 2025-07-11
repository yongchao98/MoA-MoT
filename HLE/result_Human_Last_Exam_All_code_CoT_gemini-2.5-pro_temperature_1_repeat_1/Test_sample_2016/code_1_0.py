import collections

class State:
    """
    A class to represent the state of the world using fluents.
    """
    def __init__(self):
        self.fluents = collections.defaultdict(bool)
        # Initial conditions
        self.fluents['fluent_1(object_1_type_0, object_2_type_0)'] = True
        self.fluents['fluent_1(object_2_type_0, object_1_type_0)'] = True
        self.fluents['fluent_2(object_1_type_0)'] = True

    def _apply_effects(self, add_list, del_list):
        """Applies delete effects first, then add effects."""
        for fluent in del_list:
            self.fluents[fluent] = False
        for fluent in add_list:
            self.fluents[fluent] = True

    def execute_action(self, action_name, params):
        """Executes a given action with its parameters."""
        action_methods = {
            "action_0": self.action_0,
            "action_1": self.action_1,
            "action_2": self.action_2,
            "action_3": self.action_3
        }
        if action_name in action_methods:
            action_methods[action_name](*params)
        else:
            raise ValueError(f"Unknown action: {action_name}")

    def action_0(self, p0):
        # Preconditions
        if not self.fluents[f'fluent_2({p0})']:
            raise ValueError(f"Precondition failed for action_0({p0})")
        # Effects
        self._apply_effects(add_list=[], del_list=[f'fluent_2({p0})'])

    def action_1(self, p0, p1):
        # Preconditions
        if not (not self.fluents[f'fluent_2({p0})'] and self.fluents[f'fluent_0({p0}, {p1})']):
            raise ValueError(f"Precondition failed for action_1({p0}, {p1})")
        # Effects
        add_list = [f'fluent_2({p0})', f'fluent_0({p0}, {p1})', f'fluent_1({p0}, {p1})']
        del_list = [f'fluent_2({p1})', f'fluent_1({p1}, {p0})']
        self._apply_effects(add_list, del_list)

    def action_2(self, p0, p1):
        # Preconditions
        if not self.fluents[f'fluent_2({p1})']:
            raise ValueError(f"Precondition failed for action_2({p0}, {p1})")
        # Effects
        add_list = [f'fluent_1({p0}, {p1})']
        del_list = [f'fluent_2({p1})']
        self._apply_effects(add_list, del_list)

    def action_3(self, p0, p1):
        # Preconditions
        if not (self.fluents[f'fluent_1({p1}, {p0})'] and self.fluents[f'fluent_2({p0})'] and not self.fluents[f'fluent_0({p1}, {p0})']):
            raise ValueError(f"Precondition failed for action_3({p0}, {p1})")
        # Effects
        add_list = [f'fluent_2({p0})', f'fluent_0({p1}, {p0})']
        del_list = [f'fluent_0({p0}, {p1})', f'fluent_1({p1}, {p0})']
        self._apply_effects(add_list, del_list)

    def get_final_goal(self):
        """Returns the final goal as a formatted string."""
        true_fluents = sorted([f for f, v in self.fluents.items() if v])
        return "&".join(true_fluents)

def solve():
    """
    Simulates the plan and prints the final goal.
    """
    # Define objects
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Define the plan
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

    # Initialize and run the simulation
    current_state = State()
    for action_name, params in plan:
        current_state.execute_action(action_name, params)

    # Print the final goal
    print(current_state.get_final_goal())

solve()