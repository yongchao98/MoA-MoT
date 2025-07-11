import textwrap

def explain_covalency():
    """
    Explains the relationship between orbital overlap and covalency
    and determines which compound is more covalent based on the provided data.
    """
    
    # Define the core chemical principle
    principle = "The degree of covalency in a bond is directly related to the extent of orbital overlap between the bonding atoms. Greater overlap leads to more effective electron sharing and, consequently, stronger covalent character."

    # State the observation given in the problem
    observation_CeF6 = "CeF6(2-) has greater 4f orbital overlap with F 2p orbitals."
    observation_CeCl6 = "CeCl6(2-) has lesser 4f orbital overlap with Cl 3p orbitals."

    # Apply the principle to the observations to draw a conclusion
    conclusion = "Given that CeF6(2-) has greater orbital overlap than CeCl6(2-), it will display stronger covalency."

    # Print the explanation in a clear, formatted way
    print("Chemical Principle:")
    print(textwrap.fill(principle, width=80))
    print("\nApplying the principle to the given compounds:")
    print(f"- {observation_CeF6}")
    print(f"- {observation_CeCl6}")
    print("\nConclusion:")
    print(textwrap.fill(conclusion, width=80))

explain_covalency()