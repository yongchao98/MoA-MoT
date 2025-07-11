import textwrap

def solve_welding_flow_problem():
    """
    This program analyzes the fluid dynamics in a weld pool to determine
    the dominant force mechanism based on the observed flow pattern.
    """

    # 1. Define the problem parameters from the user's query
    material = "304 stainless steel"
    observation = "outer portions of the weld pool flow inwards"
    answer_choices = {
        'A': 'Marangoni Force',
        'B': 'Arc drag Force',
        'C': 'Arc Pressure Force',
        'D': 'Lorentz (electromagnetic) Force',
        'E': 'Buoyancy Force'
    }

    # 2. Explain the core concept: The Marangoni Effect
    print("--- Analysis of Weld Pool Fluid Flow ---")
    print(f"Observation: In the weld pool of {material}, the {observation}.\n")

    explanation_marangoni = """
    The dominant mechanism for surface flow in a weld pool is often the Marangoni force, which is caused by a gradient in surface tension. Surface tension (gamma, γ) changes with temperature (T). This relationship is described by the temperature coefficient of surface tension (dγ/dT).
    """
    print(textwrap.dedent(explanation_marangoni))

    # 3. Explain the two cases of the Marangoni effect
    explanation_case1 = """
    Case 1: Negative dγ/dT (e.g., pure metals)
    - The hottest part (center) has the LOWEST surface tension.
    - The cooler part (edge) has the HIGHEST surface tension.
    - Fluid flows from low to high tension, so the flow is OUTWARD.
    - Result: A wide, shallow weld pool.
    """
    print(textwrap.dedent(explanation_case1))

    explanation_case2 = """
    Case 2: Positive dγ/dT (e.g., alloys with surfactants like sulfur)
    - The hottest part (center) has the HIGHEST surface tension.
    - The cooler part (edge) has the LOWEST surface tension.
    - Fluid flows from low to high tension, so the flow is INWARD.
    - Result: A narrow, deep weld pool.
    """
    print(textwrap.dedent(explanation_case2))

    # 4. Connect the concept to the specific problem
    print("--- Applying to the Specific Case ---")
    final_reasoning = f"""
    The material is {material}, which contains surface-active elements (surfactants) like sulfur. The presence of these elements causes the temperature coefficient of surface tension (dγ/dT) to become positive.

    Therefore, the observed INWARD flow is a direct result of the Marangoni force pulling the surface fluid from the cooler edges (lower surface tension) to the hot center (higher surface tension).
    """
    print(textwrap.dedent(final_reasoning))

    # 5. Conclude the answer
    final_answer_key = 'A'
    final_answer_text = answer_choices[final_answer_key]
    print(f"Conclusion: The dominant mechanism is the {final_answer_text}.")
    print(f"Final Answer Choice: {final_answer_key}")

solve_welding_flow_problem()