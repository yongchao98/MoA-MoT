import sys

def analyze_maxwells_demon():
    """
    Analyzes the 'Maxwell's Demon' apparatus to identify the key experimental parameter.
    """
    print("Analyzing the thought experiment of 'Maxwell's Demon'...")
    print("-" * 50)
    print("1. The Goal: To trap all gas molecules, initially in two compartments, into one.")
    print("   This process involves a decrease in the gas's entropy (it becomes more ordered).")
    print("\n2. The Second Law of Thermodynamics:")
    print("   This law states that the entropy of an isolated system at thermal equilibrium cannot decrease.")
    print("   Therefore, for the proposed sorting to work, the system must not be a simple isolated system in equilibrium.")
    print("\n3. The 'Demon' Mechanism (A One-Way Door as a Brownian Ratchet):")
    print("   The one-way door is intended to let molecules pass from compartment A to B, but not from B to A.")
    print("\n4. The Critical Flaw at Uniform Temperature (Feynman-Smoluchowski Ratchet):")
    print("   If the door (the ratchet) is at the SAME TEMPERATURE as the gas, the door itself experiences Brownian motion.")
    print("   This random thermal jiggling of the door will cause it to randomly open and allow molecules to pass in the 'wrong' direction.")
    print("   At a single, uniform temperature, the system reaches equilibrium, and no net sorting occurs.")
    print("\n5. The Required Condition for Sorting:")
    print("   To break this symmetry and make the door a true one-way gate, it must NOT be in thermal equilibrium with the gas.")
    print("   Specifically, the door mechanism must be kept at a much LOWER TEMPERATURE than the gas.")
    print("   - A cold door does not jiggle randomly, so it can perform its sorting function effectively.")
    print("   - This requires an external energy input to keep the door cold, meaning the entire system (gas + door + cooling system) does obey the Second Law.")
    print("-" * 50)
    print("\nConclusion:")
    print("The essential experimental parameter that enables the sorting process by creating a non-equilibrium condition is TEMPERATURE.")

analyze_maxwells_demon()