import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying physical constraints.
    
    The function checks for:
    1. Correct scaling with the symmetry-breaking VEVs.
    2. Completeness of the particle spectrum in the loop calculation.
    3. Correct signs for boson (+) and fermion (-) contributions.
    """
    
    # The LLM's answer to be checked
    llm_answer = 'A'

    # The options from the question, using raw strings for regex compatibility
    options = {
        'A': r"M_{h_{2}}^{2}=\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}-\alpha_{4}M_{t}^{4}+\alpha_{5}M_{H^{\pm}}^{4}+\alpha_{6}M_{H^{0}}^{4}+\alpha_{7}M_{A^{0}}^{4}-\alpha_{8}\sum M_{N_{i}}^{4}\right\}",
        'B': r"M_{h_{2}}^{2}=\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}-\alpha_{4}M_{t}^{4}+\alpha_{5}M_{H^{\pm}}^{4}+\alpha_{6}M_{H^{0}}^{4}-\alpha_{7}\sum M_{N_{i}}^{4}\right\}",
        'C': r"M_{h_{2}}^{2}=\frac{\left(x^{2}+v^{2}\right)}{8\pi^{2}}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}-\alpha_{4}M_{t}^{4}+\alpha_{5}M_{H^{\pm}}^{4}+\alpha_{6}M_{H^{0}}^{4}+\alpha_{7}M_{A^{0}}^{4}-\alpha_{8}\sum M_{N_{i}}^{4}\right\}",
        'D': r"M_{h_{2}}^{2}=\frac{1}{8\pi^{2}\left(x^{2}+v^{2}\right)}\left\{ \alpha_{1}M_{h_{1}}^{4}+\alpha_{2}M_{W}^{4}+\alpha_{3}M_{Z}^{4}+\alpha_{4}M_{H^{\pm}}^{4}+\alpha_{5}M_{H^{0}}^{4}+\alpha_{6}M_{A^{0}}^{4}-\alpha_{7}\sum M_{N_{i}}^{4}\right\}"
    }

    # Define the expected particle content and their types
    bosons = ["M_{h_{1}}^{4}", "M_{W}^{4}", "M_{Z}^{4}", "M_{H^{\pm}}^{4}", "M_{H^{0}}^{4}", "M_{A^{0}}^{4}"]
    fermions = ["M_{t}^{4}", r"\sum M_{N_{i}}^{4}"]
    all_particles = bosons + fermions

    correct_option_key = None
    error_log = {}

    for key, formula in options.items():
        errors = []
        
        # Constraint 1: Check scaling factor (must be in denominator)
        if r"\frac{\left(x^{2}+v^{2}\right)" in formula:
            errors.append("Incorrect scaling: The factor (x^2+v^2) must be in the denominator.")
        
        # Constraint 2: Check for completeness of particle content
        missing_particles = [p for p in all_particles if p not in formula]
        if missing_particles:
            # Clean up names for printing
            clean_names = [p.replace('^{4}', '').replace(r'\sum ', '') for p in missing_particles]
            errors.append(f"Incomplete particle spectrum: Missing contributions from {', '.join(clean_names)}.")

        # Constraint 3: Check signs of contributions
        # Fermions must have a negative sign
        for fermion in fermions:
            pattern = r"([+-])\s*\\alpha_{\d+}\s*" + re.escape(fermion)
            match = re.search(pattern, formula)
            if not match or match.group(1) != '-':
                clean_name = fermion.replace('^{4}', '').replace(r'\sum ', '')
                errors.append(f"Incorrect sign for fermion '{clean_name}'. It must be negative.")
        
        # Bosons must have a positive sign (i.e., not be preceded by a minus)
        for boson in bosons:
            pattern = r"([+-])\s*\\alpha_{\d+}\s*" + re.escape(boson)
            match = re.search(pattern, formula)
            if match and match.group(1) == '-':
                clean_name = boson.replace('^{4}', '')
                errors.append(f"Incorrect sign for boson '{clean_name}'. It must be positive.")

        if not errors:
            correct_option_key = key
        else:
            error_log[key] = errors

    if llm_answer == correct_option_key:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect."
        if llm_answer in error_log:
            reason += f" The following constraints are not satisfied by option {llm_answer}:\n"
            for err in error_log[llm_answer]:
                reason += f"- {err}\n"
        else:
             reason += f" The correct answer is '{correct_option_key}' because it is the only option that satisfies all physical constraints."
        return reason

# Execute the check and print the result
result = check_correctness()
print(result)