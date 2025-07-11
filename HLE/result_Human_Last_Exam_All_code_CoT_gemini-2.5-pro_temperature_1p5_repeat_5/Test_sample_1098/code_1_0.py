import textwrap

def explain_maxwell_demon_ratchet():
    """
    Explains why temperature is the key parameter for a Brownian ratchet
    to concentrate gas on one side of a two-compartment chamber.
    """

    # The core of the problem is the Second Law of Thermodynamics.
    # Moving all gas to one side is a decrease in entropy (disorder).
    # entropy_final << entropy_initial
    # This cannot happen spontaneously in an isolated system at equilibrium.
    explanation = """
    The problem describes a physical version of the Maxwell's demon thought experiment using a Brownian ratchet. The goal is to collect all gas particles into one of two compartments.

    1.  **Entropy and the Second Law of Thermodynamics:** Spreading gas across two compartments is a high-disorder (high-entropy) state. Confining all the gas to one compartment is a low-disorder (low-entropy) state. The Second Law of Thermodynamics states that the total entropy of an isolated system cannot decrease. Therefore, this process cannot happen spontaneously if the entire system is at a single, uniform temperature (in thermal equilibrium).

    2.  **The Brownian Ratchet:** A one-way door acting as a Brownian ratchet might seem to work. However, if the door is at the same temperature as the gas, the door itself will experience random thermal jiggling. This random motion will cause the 'one-way' mechanism to fail, allowing particles to move in both directions at an equal rate. There will be no net flow.

    3.  **The Required Condition:** To make the ratchet work and pump the gas to one side, there must be a temperature difference. This temperature difference drives the system out of equilibrium. For example, if the ratchet mechanism is kept colder than the gas, it can rectify the random motion of the gas particles. A temperature difference provides the necessary energy gradient to perform the work of compressing the gas, which is consistent with the laws of thermodynamics.

    4.  **Conclusion:** Of the choices provided, only Temperature creates the necessary non-equilibrium condition for the device to function. Therefore, a temperature difference is the required experimental parameter.
    """

    print(textwrap.dedent(explanation).strip())

explain_maxwell_demon_ratchet()