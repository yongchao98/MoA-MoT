import sys

def analyze_evolutionary_challenge():
    """
    This function models the challenge to genetic drift predictions
    posed by adaptive evolution.
    """
    # Predictive models of genetic drift are based on the neutral theory,
    # where the nonsynonymous substitution rate (Ka) is expected to be
    # roughly equal to the synonymous substitution rate (Ks).
    #
    # Ks = Rate of synonymous (silent) substitutions, our proxy for the neutral/drift rate.
    # Ka = Rate of nonsynonymous (amino-acid changing) substitutions.

    # We will simulate a scenario of adaptive evolution.
    # In this scenario, selection favors new protein variants, causing
    # nonsynonymous changes to accumulate faster than neutral mutations.
    ka_observed = 0.12  # Observed nonsynonymous rate is high
    ks_observed = 0.04  # Observed synonymous (neutral) rate

    print("Analyzing a scenario that challenges genetic drift models...")
    print(f"Observed nonsynonymous rate (Ka): {ka_observed}")
    print(f"Observed synonymous rate (Ks): {ks_observed}")
    print("-" * 30)

    # Under adaptive evolution, the key relationship is Ka > Ks.
    # The models of genetic drift are challenged because they predict Ka ≈ Ks.
    # Let's check this condition.
    
    if ka_observed > ks_observed:
        ratio = ka_observed / ks_observed
        print(f"The ratio Ka/Ks is {ratio:.2f}.")
        print("This value is significantly greater than 1.")
        print("\nConclusion: The increased variability in nonsynonymous sites (high Ka)")
        print("driven by adaptive evolution provides a result that outweighs and contradicts")
        print("the predictions from models based purely on random genetic drift.")

        # Per the instructions, printing each number in the final conceptual equation.
        # The conceptual equation is a test for positive selection: Ka / Ks > 1
        print("\nBreaking down the challenge equation (Ka / Ks > 1):")
        print(f"Component 'Ka': {ka_observed}")
        print(f"Component 'Ks': {ks_observed}")
        print(f"Comparison threshold: 1")
    else:
        print("The observed rates are consistent with neutral evolution or purifying selection.")
        print("This does not strongly challenge genetic drift models.")

analyze_evolutionary_challenge()