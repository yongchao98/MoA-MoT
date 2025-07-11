import sys

# Python's recursion limit can be a constraint for deep justification chains.
# We can increase it for demonstration, though our example uses a direct loop.
sys.setrecursionlimit(2000)

# This dictionary represents the objective, external truth of the world.
# The agent does NOT have direct access to this; it's the "God's-eye view".
world_state = {
    "The sun is a star": True,
    "The sky is blue due to Rayleigh scattering": True,
    "My keys are on the kitchen table": False, # In reality, they fell behind the couch.
    "I saw myself put the keys on the table": True, # This is a true memory...
    "The table is in the kitchen": True, # ...based on true premises...
    "Therefore, my keys are on the kitchen table": False # ...but the conclusion is false.
}

# This dictionary represents the agent's internal web of beliefs and their justifications.
# A belief is justified if the propositions in its justification set are known.
belief_system = {
    "The sun is a star": {"The sky is blue due to Rayleigh scattering"},
    # Problem 1: This is a circular justification. To know A, you must know B, and to know B, you must know A.
    "The sky is blue due to Rayleigh scattering": {"The sun is a star"},

    # Problem 2: A well-justified belief that is false.
    "My keys are on the kitchen table": {
        "I saw myself put the keys on the table",
        "The table is in the kitchen"
    }
}

def has_knowledge(proposition, justification_path=None):
    """
    Checks if a proposition constitutes 'knowledge' according to the JTB definition.
    This function has access to the external `world_state` to check for Truth.
    
    Args:
        proposition (str): The belief being checked.
        justification_path (set): Used internally to track the chain of justifications to detect loops.
    """
    if justification_path is None:
        justification_path = set()

    print(f"\n--- Checking knowledge of: '{proposition}' ---")

    # 1. Belief (B): We assume if it's in the belief_system, the agent believes it.
    if proposition not in belief_system:
        print(f"CONCLUSION: Agent does NOT have knowledge. Proposition not in belief system.")
        return False
    print(f"B (Belief): The agent holds the belief '{proposition}'.")

    # 2. Truth (T): Check against the external `world_state`. The agent cannot do this.
    is_true = world_state.get(proposition, False)
    print(f"T (Truth): Is '{proposition}' objectively true? ==> {is_true}")
    if not is_true:
        # This demonstrates Problem 2: The Inaccessibility of Truth.
        # The agent's justification might be flawless, but if the belief is false, it isn't knowledge.
        print(f"CONCLUSION: Agent does NOT have knowledge of '{proposition}'. The belief is FALSE.")
        return False

    # 3. Justification (J): Check the agent's reasons for the belief.
    print(f"J (Justification): Checking justification for '{proposition}'...")

    # This check demonstrates Problem 1: The Regress/Circularity Problem.
    if proposition in justification_path:
        print(f"JUSTIFICATION FAILURE: Circular reasoning detected! We are already trying to justify '{proposition}' in this chain.")
        print(f"CONCLUSION: Agent does NOT have knowledge of '{proposition}'. The justification is circular.")
        return False
    justification_path.add(proposition)

    justifiers = belief_system.get(proposition)
    if not justifiers: # It's a foundational belief
        print(f"Justification: '{proposition}' is a foundational belief. As it's true, we'll accept it as justified.")
        # In a real system, how to validate foundational beliefs is the core of the problem.
        justification_path.remove(proposition)
        return True

    print(f"The justification for '{proposition}' requires knowledge of: {justifiers}")
    
    # Recursively check if the justifying beliefs are themselves knowledge.
    for justifier in justifiers:
        if not has_knowledge(justifier, justification_path):
            print(f"JUSTIFICATION FAILURE: The justifying belief '{justifier}' is not known.")
            print(f"CONCLUSION: Agent does NOT have knowledge of '{proposition}'.")
            justification_path.remove(proposition)
            return False

    # If all checks pass...
    justification_path.remove(proposition)
    print(f"CONCLUSION: Agent HAS KNOWLEDGE of '{proposition}'. The belief is Justified and True.")
    return True

# --- Main Demonstration ---
print("="*60)
print("Demonstrating Problem 1: The Regress/Circularity Problem")
print("="*60)
print("We check 'The sun is a star', which has a circular justification.")
has_knowledge("The sun is a star")

print("\n\n" + "="*60)
print("Demonstrating Problem 2: The Inaccessibility of Truth")
print("="*60)
print("We check 'My keys are on the kitchen table'. The agent has good reasons")
print("to believe this, but the belief is objectively false.")
# We need to add the justifying beliefs to the system as foundational.
belief_system["I saw myself put the keys on the table"] = set()
belief_system["The table is in the kitchen"] = set()
has_knowledge("My keys are on the kitchen table")