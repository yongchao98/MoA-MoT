def check_answer_correctness():
    """
    Checks the correctness of the proposed answer for the given chemistry question.

    The question is: Identify the product when (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
    is reacted with Me2CuLi.

    The proposed final answer from the LLM is:
    C) (1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol
    """

    # Define the constraints based on the problem description and chemical principles.
    # Constraint 1: Regioselectivity - Attack at the less hindered carbon (C6).
    # This leads to a 1,2,4,5-tetramethylcyclohexan-1-ol skeleton.
    correct_skeleton = "1,2,4,5-tetramethylcyclohexan-1-ol"

    # Constraint 2: Stereoselectivity - Inversion of configuration at the attacked carbon.
    # Starting C6 configuration is (S). Inversion means the product C2 configuration should be (R).
    expected_config_at_C2 = "2R"

    # Constraint 3: Retention of configuration at other centers.
    # Starting C1 is (R), C3 is (R), C4 is (R).
    # Product should have (1R), (5R from C3), (4R).
    expected_retained_configs = ["1R", "4R", "5R"]

    # The answer to check.
    llm_answer = "C) (1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"

    # --- Perform Checks ---

    # 1. Check skeleton/regioselectivity
    if correct_skeleton not in llm_answer:
        return (f"Incorrect. The regioselectivity is wrong. The product skeleton should be "
                f"'{correct_skeleton}' based on attack at the less hindered carbon C6.")

    # 2. Check stereochemistry
    # Extract stereochemical descriptors from the answer string
    try:
        configs_str = llm_answer.split(')')[0].split('(')[1]
        provided_configs = [c.strip() for c in configs_str.split(',')]
    except IndexError:
        return "Incorrect. The answer format is invalid and stereochemistry cannot be parsed."

    # Check retention
    for config in expected_retained_configs:
        if config not in provided_configs:
            return (f"Incorrect. The configuration at an unreacted stereocenter is wrong. "
                    f"Expected to find '{config}' but it was not present in the answer '{llm_answer}'.")

    # Check inversion
    provided_config_at_C2 = next((c for c in provided_configs if c.startswith('2')), None)
    if not provided_config_at_C2:
         return (f"Incorrect. The stereochemistry at C2 is missing in the answer '{llm_answer}'.")

    if provided_config_at_C2 != expected_config_at_C2:
        return (f"Incorrect. The answer violates the stereoselectivity constraint. "
                f"The reaction involves 'inversion of configuration' at the attacked carbon (C6), which has an (S) configuration. "
                f"This should result in an (R) configuration at the new C2 position. "
                f"The answer proposes an (S) configuration ({provided_config_at_C2}), which is incorrect.")

    # If all checks pass
    return "Correct"

# Run the check and print the result
result = check_answer_correctness()
if result != "Correct":
    print(result)
else:
    print(result)