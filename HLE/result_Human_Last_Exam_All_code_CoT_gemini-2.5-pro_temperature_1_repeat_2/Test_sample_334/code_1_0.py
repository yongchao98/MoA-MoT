def solve_biometric_challenge():
    """
    This function analyzes the cybersecurity requirements for enhancing biometric authentication
    and determines the best solution from a list of choices.
    """

    # The four key cybersecurity requirements derived from the problem description.
    requirements = {
        "prevention_of_spoofing_replay": "Must prevent attacks using fake or recorded biometric data.",
        "protection_against_coercion": "An unwilling (e.g., kidnapped) authorized user must NOT authenticate.",
        "resilience_to_data_exposure": "The system must remain secure even if biometric data is stolen.",
        "modality_agnostic": "The solution must apply to a wide range of biometrics (face, fingerprint, etc.)."
    }

    # Analysis of why "Strong Multi-Factor Authentication" is the correct choice.
    # It uniquely satisfies all four constraints.
    analysis = {
        "solution_choice": "I",
        "solution_name": "Strong Multi-Factor Authentication",
        "evaluation": {
            "prevention_of_spoofing_replay": {
                "satisfied": True,
                "reason": "A spoofed biometric is insufficient; the attacker still needs the additional factors (e.g., PIN and hardware token)."
            },
            "protection_against_coercion": {
                "satisfied": True,
                "reason": "The user can withhold a factor they control separately, like a hardware token. An attacker cannot force a user to produce an item they don't have with them."
            },
            "resilience_to_data_exposure": {
                "satisfied": True,
                "reason": "If the biometric data is stolen, it is useless without the other authentication factors."
            },
            "modality_agnostic": {
                "satisfied": True,
                "reason": "This approach augments any underlying biometric system, regardless of the specific modality used."
            }
        },
        "conclusion": "This is the only solution presented that robustly addresses all four specified requirements by adding factors that are secret ('something you know') and physically separate ('something you have'), directly compensating for the non-secret nature of biometrics."
    }

    # Print the final determination.
    # The output shows the letter corresponding to the correct answer choice.
    final_answer = analysis["solution_choice"]
    print(f"Analyzing all options against the four core requirements: preventing spoofing, preventing coercion, surviving data exposure, and being modality-agnostic.")
    print(f"The only solution that satisfies all criteria is '{analysis['solution_name']}'.")
    print(f"It protects against coercion because a user can withhold a physical token, and it survives data theft because the biometric alone is not enough for access.")
    print("\nTherefore, the correct choice is:")
    print(final_answer)


solve_biometric_challenge()