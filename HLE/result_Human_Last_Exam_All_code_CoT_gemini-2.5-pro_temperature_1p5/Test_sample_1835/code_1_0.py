# 1. Define a "world" of objects and their properties.
# These represent the things we can think about.
world = {
    'apple': {'color': 'red', 'shape': 'round', 'is_heavy': False},
    'brick': {'color': 'red', 'shape': 'rectangular', 'is_heavy': True},
    'sky': {'color': 'blue', 'shape': 'vast', 'is_heavy': False},
}

# 2. Define predicates (concepts like 'F', 'G') as Python functions.
# Understanding a proposition like "the apple is red" means grasping a function like `is_red`.
def is_red(item_name):
    """This function represents the predicate concept 'is_red'."""
    return world.get(item_name, {}).get('color') == 'red'

def is_heavy(item_name):
    """This function represents the predicate concept 'is_heavy'."""
    return world.get(item_name, {}).get('is_heavy', False)

# 3. Define the concepts we "understand": subjects and predicates.
subjects = list(world.keys())
predicates = {
    "is_red": is_red,
    "is_heavy": is_heavy
}

print("--- Premise: We understand individual propositions ---")
# Example: Understanding 'the apple is red' (Fa)
subject_a = 'apple'
print(f"We can think '{is_red.__name__}(\"{subject_a}\")'. This yields: {is_red(subject_a)}")
# Example: Understanding 'the brick is heavy' (Gb)
subject_b = 'brick'
print(f"We can think '{is_heavy.__name__}(\"{subject_b}\")'. This yields: {is_heavy(subject_b)}\n")


print("--- Demonstrating the Generality Constraint (Recombinability) ---")
print("Because we grasp the concepts separately, we can form new thoughts:")
for p_name, p_func in predicates.items():
    for s_name in subjects:
        print(f"Formed thought: '{p_name}(\"{s_name}\")' -> {p_func(s_name)}")
print("-" * 25 + "\n")


print("--- Extending to Universal Quantification ---")
print("We also understand universal quantification (∀x), modeled here as the 'for_all' function.")
print("We can combine a known predicate (e.g., 'is_red') with 'for_all' to form and evaluate '∀x, is_red(x)'.")

def for_all(predicate_func, domain):
    """
    This function models the concept of universal quantification (∀x).
    It checks if a predicate holds true for all items in a domain.
    """
    print(f"\nEVALUATING: 'For all x in {domain}, does {predicate_func.__name__}(x) hold?'")
    print("The components of this proposition are:")
    
    final_result = True
    for item in domain:
        # Here we evaluate each component of the final "equation"
        component_result = predicate_func(item)
        print(f"  {predicate_func.__name__}('{item}') = {component_result}")
        if not component_result:
            final_result = False
            
    print(f"THE FINAL RESULT IS: {final_result}")
    return final_result

# Now, we combine the predicate concept 'is_red' with our understanding of universal quantification.
for_all(is_red, subjects)