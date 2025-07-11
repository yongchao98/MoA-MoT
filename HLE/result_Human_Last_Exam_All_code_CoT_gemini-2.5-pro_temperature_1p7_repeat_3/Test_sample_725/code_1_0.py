def analyze_covalency():
    """
    Analyzes the relative covalency of CeF6(2-) and CeCl6(2-) based on orbital overlap.
    """

    # 1. Define the key information from the problem description.
    compound_1 = "CeF6²⁻"
    compound_2 = "CeCl6²⁻"
    observation = "The 4f-2p orbital overlap in CeF6²⁻ is greater than the 4f-3p overlap in CeCl6²⁻."

    # 2. State the fundamental chemical principle connecting orbital overlap and covalency.
    principle = "In chemical bonding, covalency arises from the sharing of electrons via the overlap of atomic orbitals. A greater degree of orbital overlap signifies more effective electron sharing and therefore corresponds to a stronger covalent character."

    # 3. Apply the principle to the provided observation to draw a conclusion.
    conclusion = f"Given that {compound_1} has greater orbital overlap than {compound_2}, it logically follows from the chemical principle that {compound_1} will exhibit stronger covalency."

    # 4. Print the step-by-step reasoning.
    print("Step 1: Identify the core chemical principle.")
    print(f"Principle: {principle}\n")

    print("Step 2: Analyze the given information.")
    print(f"Observation: {observation}\n")

    print("Step 3: Draw a conclusion by applying the principle to the observation.")
    print(f"Conclusion: {conclusion}\n")

    print("Final Answer:")
    print("Based on its greater orbital overlap, CeF6²⁻ would display stronger covalency compared to CeCl6²⁻.")


if __name__ == "__main__":
    analyze_covalency()