def solve_geology_question():
    """
    Analyzes the geological transition from Rock A to Rock B and determines the plausible causes.
    """
    # Step 1: Define the geological interpretation.
    # The transition from a dark, fine-grained rock (A, deep water) to a lighter,
    # thicker-bedded rock (B, shallower water) signifies a marine regression (fall in sea level).
    interpretation = "marine regression"
    print(f"Geological Interpretation: The transition from Rock A to Rock B is a 'shallowing-upward' sequence, indicating a {interpretation}.")
    print("-" * 30)

    # Step 2: Define the options and their effects on sea level.
    options = {
        'i': {'event': 'Onset of glaciation', 'effect': 'sea level fall', 'result': 'marine regression'},
        'ii': {'event': 'Coastal subsidence', 'effect': 'relative sea level rise', 'result': 'marine transgression'},
        'iii': {'event': 'Global warming', 'effect': 'sea level rise', 'result': 'marine transgression'},
        'iv': {'event': 'Decreased Mid-Ocean Ridge Activity', 'effect': 'sea level fall', 'result': 'marine regression'}
    }

    # Step 3: Identify the plausible options by matching the result with the interpretation.
    plausible_options = []
    print("Evaluating the options:")
    for key, value in options.items():
        is_plausible = value['result'] == interpretation
        plausibility_text = "Plausible" if is_plausible else "Not Plausible"
        print(f"  {key}. {value['event']}: Causes a {value['effect']}, leading to a {value['result']}. --> {plausibility_text}")
        if is_plausible:
            plausible_options.append(key)

    # Step 4: Format the final answer as a comma-separated string.
    final_answer = ",".join(plausible_options)
    print("-" * 30)
    print(f"The plausible explanations are those that cause a marine regression.")
    print(f"Final Answer (comma-separated list): {final_answer}")

solve_geology_question()