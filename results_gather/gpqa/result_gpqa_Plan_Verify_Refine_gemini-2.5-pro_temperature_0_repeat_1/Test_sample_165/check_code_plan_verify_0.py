import re

def check_physics_answer():
    """
    Checks the correctness of the LLM's answer by verifying its reasoning
    based on the Gildener-Weinberg formalism.
    """

    # --- Define Physics Rules and Model Particles ---

    # Rule 1: The correct prefactor for the mass-squared formula.
    CORRECT_PREFACTOR = "1 / (8*pi^2 * (x^2+v^2))"
    
    # Rule 2: Particle types and the expected sign of their M^4 contribution.
    # Bosons contribute positively (+), Fermions contribute negatively (-).
    PARTICLE_PROPERTIES = {
        "h1":   {"type": "boson", "expected_sign": "+"},
        "W":    {"type": "boson", "expected_sign": "+"},
        "Z":    {"type": "boson", "expected_sign": "+"},
        "t":    {"type": "fermion", "expected_sign": "-"},
        "H_pm": {"type": "boson", "expected_sign": "+"}, # H^\pm
        "H0":   {"type": "boson", "expected_sign": "+"},
        "A0":   {"type": "boson", "expected_sign": "+"},
        "N":    {"type": "fermion", "expected_sign": "-"}, # N_i
    }

    # Rule 3: The complete set of particles expected to contribute.
    EXPECTED_PARTICLES = set(PARTICLE_PROPERTIES.keys())

    # --- Problem Data ---

    # The options provided in the question.
    # Using raw strings and simplified notation for parsing.
    options_str = {
        "A": r"M_{h_{2}}^{2}=\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}+\alpha_{4}M_{H^{\pm}}^{4}+\alpha_{5}M_{H^{0}}^{4}+\alpha_{6}M_{A^{0}}^{4}-\alpha_{7}\sum M_{N_{i}}^{4}\right\}",
        "B": r"M_{h_{2}}^{2}=\frac{\left(x^{2}+v^{2}\right)}{8\pi^{2}}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}-\alpha_{4}M_{t}^{4}+\alpha_{5}M_{H^{\pm}}^{4}+\alpha_{6}M_{H^{0}}^{4}+\alpha_{7}M_{A^{0}}^{4}-\alpha_{8}\sum M_{N_{i}}^{4}\right\}",
        "C": r"M_{h_{2}}^{2}=\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}-\alpha_{4}M_{t}^{4}+\alpha_{5}M_{H^{\pm}}^{4}+\alpha_{6}M_{H^{0}}^{4}-\alpha_{7}\sum M_{N_{i}}^{4}\right\}",
        "D": r"M_{h_{2}}^{2}=\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}-\alpha_{4}M_{t}^{4}+\alpha_{5}M_{H^{\pm}}^{4}+\alpha_{6}M_{H^{0}}^{4}+\alpha_{7}M_{A^{0}}^{4}-\alpha_{8}\sum M_{N_{i}}^{4}\right\}"
    }
    
    # The LLM's final answer choice.
    llm_choice = "D"

    # --- Parsing and Verification Logic ---

    def parse_formula(formula_str):
        """Parses a formula string into a structured dictionary."""
        parsed = {"prefactor": "unknown", "terms": {}}
        
        # 1. Parse Prefactor
        if r"\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}" in formula_str:
            parsed["prefactor"] = CORRECT_PREFACTOR
        elif r"\frac{\left(x^{2}+v^{2}\right)}{8\pi^{2}}" in formula_str:
            parsed["prefactor"] = "(x^2+v^2) / (8*pi^2)"

        # 2. Parse Terms
        term_map = {
            'M_{h_{1}}^{4}': 'h1', 'M_{W}^{4}': 'W', 'M_{Z}^{4}': 'Z',
            'M_{t}^{4}': 't', r'M_{H^{\pm}}^{4}': 'H_pm', 'M_{H^{0}}^{4}': 'H0',
            'M_{A^{0}}^{4}': 'A0', r'M_{N_{i}}^{4}': 'N'
        }
        
        content_match = re.search(r'\\left\{ (.*) \\right\}', formula_str)
        if not content_match: return parsed
        content = content_match.group(1)
        
        # Add implicit '+' to the first term
        if not content.strip().startswith(('+', '-')):
            content = '+' + content

        # Find all terms and their preceding signs
        terms = re.findall(r'([+-])\s*(?:\\alpha_{\d+})?\s*(?:\\sum)?\s*([\w{\\^}]+)', content)
        
        for sign, term_str in terms:
            for latex, key in term_map.items():
                if latex == term_str:
                    parsed["terms"][key] = sign
                    break
        return parsed

    # --- Execute Checks ---

    # The LLM claims D is correct and A, B, C are wrong for specific reasons.
    # We verify these claims.
    
    parsed_options = {key: parse_formula(val) for key, val in options_str.items()}

    # Verify reason for rejecting B (prefactor)
    if parsed_options["B"]["prefactor"] == CORRECT_PREFACTOR:
        return "Incorrect. The reasoning that Option B has the wrong prefactor is flawed, as the parser found the correct prefactor."
    
    # Verify reason for rejecting A (missing top quark)
    if "t" in parsed_options["A"]["terms"]:
        return "Incorrect. The reasoning that Option A is wrong because it's missing the top quark contribution is flawed. The term was found."
    if "t" not in parsed_options["D"]["terms"]: # Sanity check
        return "Incorrect. The supposedly correct Option D is missing the top quark."

    # Verify reason for rejecting C (missing pseudoscalar A0)
    if "A0" in parsed_options["C"]["terms"]:
        return "Incorrect. The reasoning that Option C is wrong because it's missing the A0 contribution is flawed. The term was found."
    if "A0" not in parsed_options["D"]["terms"]: # Sanity check
        return "Incorrect. The supposedly correct Option D is missing the pseudoscalar A0."

    # Now, fully check the chosen answer, D.
    chosen_data = parsed_options[llm_choice]

    # Check 1: Prefactor
    if chosen_data["prefactor"] != CORRECT_PREFACTOR:
        return f"Incorrect. The final answer {llm_choice} is wrong because its prefactor is '{chosen_data['prefactor']}' but should be '{CORRECT_PREFACTOR}'."

    # Check 2: Completeness
    present_particles = set(chosen_data["terms"].keys())
    if present_particles != EXPECTED_PARTICLES:
        missing = EXPECTED_PARTICLES - present_particles
        extra = present_particles - EXPECTED_PARTICLES
        error = f"Incorrect. The final answer {llm_choice} is wrong because it has an incomplete set of particles."
        if missing: error += f" It is missing: {sorted(list(missing))}."
        if extra: error += f" It incorrectly includes: {sorted(list(extra))}."
        return error

    # Check 3: Signs
    for particle, sign in chosen_data["terms"].items():
        props = PARTICLE_PROPERTIES[particle]
        if sign != props["expected_sign"]:
            return (f"Incorrect. The final answer {llm_choice} is wrong. "
                    f"The particle '{particle}' is a {props['type']} and should have a '{props['expected_sign']}' sign, "
                    f"but it has a '{sign}' sign in the formula.")

    # If all checks on the reasoning and the final answer pass
    return "Correct"

# Run the check and print the result
result = check_physics_answer()
print(result)