import collections

def check_maxwell_equations_answer():
    """
    Checks the correctness of an answer about which of Maxwell's equations
    would change if magnetic monopoles existed.
    """

    # 1. Define the core concepts of Maxwell's equations.
    # The keys are the physical descriptions used in the question's options.
    # The values are the standard equations and a boolean indicating if they
    # change with the introduction of magnetic monopoles.
    # A magnetic monopole is a source of the magnetic field (like an electric charge
    # is for the electric field), and its movement constitutes a magnetic current.
    maxwell_info = {
        # Gauss's Law for Electricity: Unchanged. Deals with electric charges.
        "divergence of the electric field": {"changes": False, "reason": "This law relates the electric field to electric charges, which are unaffected by magnetic monopoles."},
        
        # Gauss's Law for Magnetism: MUST change. The standard law (∇ ⋅ B = 0) is the
        # mathematical statement that monopoles DON'T exist.
        "divergence of the magnetic field": {"changes": True, "reason": "The standard law (∇ ⋅ B = 0) states there are no magnetic monopoles. If they existed, they would act as sources for the magnetic field, requiring a magnetic charge density term (∇ ⋅ B = ρ_m)."},
        
        # Faraday's Law of Induction: MUST change. By symmetry, a moving magnetic monopole
        # (a magnetic current) must induce a circulating electric field.
        "circulation of the electric field": {"changes": True, "reason": "The standard law (∇ × E = -∂B/∂t) would need an additional term for magnetic current (J_m), analogous to how electric current creates a magnetic field. The new law would be ∇ × E = -∂B/∂t - J_m."},
        
        # Ampère-Maxwell Law: Unchanged. Deals with electric currents and changing
        # electric fields as sources for the magnetic field.
        "circulation of the magnetic field": {"changes": False, "reason": "This law relates the magnetic field to electric currents and changing electric fields, which are not directly altered by the existence of magnetic monopoles."}
    }

    # 2. Define the options from the question.
    # We normalize terms like "flux" to "divergence" and "curl" to "circulation"
    # for consistent matching.
    options = {
        "A": {"circulation of the magnetic field", "divergence of the electric field"}, # "flux of the electric field" -> divergence
        "B": {"divergence of the magnetic field"},
        "C": {"circulation of the electric field", "divergence of the magnetic field"},
        "D": {"divergence of the magnetic field", "circulation of the magnetic field"} # "curl of the magnetic field" -> circulation
    }

    # The answer provided by the LLM
    llm_answer_choice = "C"

    # 3. Determine the ground truth from our physics model.
    # This is the set of descriptions for equations that are expected to change.
    ground_truth_changes = {desc for desc, info in maxwell_info.items() if info["changes"]}

    # 4. Get the set of concepts from the LLM's chosen answer.
    if llm_answer_choice not in options:
        return f"Invalid answer choice '{llm_answer_choice}'. The choice must be one of {list(options.keys())}."
        
    llm_answer_concepts = options[llm_answer_choice]

    # 5. Compare the LLM's answer with the ground truth.
    if llm_answer_concepts == ground_truth_changes:
        return "Correct"
    else:
        # Analyze the discrepancy to provide a clear reason.
        missing_from_answer = ground_truth_changes - llm_answer_concepts
        incorrectly_in_answer = llm_answer_concepts - ground_truth_changes

        error_messages = []
        if missing_from_answer:
            for concept in missing_from_answer:
                error_messages.append(f"The answer is incomplete. It fails to mention that the law for the '{concept}' must also change. Reason: {maxwell_info[concept]['reason']}")
        
        if incorrectly_in_answer:
            for concept in incorrectly_in_answer:
                error_messages.append(f"The answer is incorrect. It wrongly includes the law for the '{concept}' as changing. Reason: {maxwell_info[concept]['reason']}")
        
        return "\n".join(error_messages)

# Run the check
result = check_maxwell_equations_answer()
print(result)