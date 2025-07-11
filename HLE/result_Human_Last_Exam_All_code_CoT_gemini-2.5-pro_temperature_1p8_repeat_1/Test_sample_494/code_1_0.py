import textwrap

def explain_ammonia_tunneling():
    """
    Explains whether an exotic form of ammonia would exhibit tunneling.
    The exotic ammonia has hydrogens with spin-0 nuclei (bosons) instead of
    the usual spin-1/2 nuclei (fermions).
    """

    # --- Introduction ---
    print("--- The Question: Ammonia Tunneling with Exotic Hydrogens ---")
    explanation = """
    This script analyzes whether an ammonia molecule (NH3) would still exhibit its famous inversion tunneling if its three ordinary hydrogen atoms were replaced with 'exotic' hydrogen atoms.
    
    - Ordinary Hydrogen: Nucleus is a proton (spin 1/2, a fermion).
    - Exotic Hydrogen:   Nucleus has spin 0 (a boson).
    """
    print(textwrap.dedent(explanation))

    # --- Step 1: What causes tunneling in ammonia? ---
    print("\n--- Step 1: The Physical Basis of Tunneling ---")
    explanation = """
    Ammonia has a pyramid shape. The Nitrogen atom can sit on one side of the plane formed by the three hydrogens or on the other. These two states are separated by a potential energy barrier.
    
    Quantum Tunneling is the process that allows the Nitrogen atom to pass *through* this barrier, even without enough energy to go over it. The existence of tunneling depends directly on the existence of this finite potential energy barrier.
    """
    print(textwrap.dedent(explanation))

    # --- Step 2: Does changing nuclear spin affect the barrier? ---
    print("\n--- Step 2: Analyzing the Effect of Changing Nuclear Spin ---")
    explanation = """
    The potential energy barrier is created by the electrostatic forces between all the electrons and nuclei in the molecule. The key factors are the particles' charges and their geometric arrangement.
    
    In this problem, the exotic hydrogen's nucleus still has the same charge (+1) as a regular proton. Changing the nuclear spin from 1/2 to 0 has no effect on these electrostatic forces.
    
    Conclusion: The potential energy barrier remains UNCHANGED.
    """
    print(textwrap.dedent(explanation))

    # --- Step 3: What does change? ---
    print("\n--- Step 3: The Role of Quantum Statistics (Pauli Principle) ---")
    explanation = """
    While the barrier doesn't change, the quantum statistics of the system do.
    
    - For ordinary NH3 (spin-1/2, fermion hydrogens), the total wavefunction must be ANTISYMMETRIC when two hydrogens are exchanged.
    - For exotic NH3 (spin-0, boson hydrogens), the total wavefunction must be SYMMETRIC when two hydrogens are exchanged.
    
    This rule affects which rotational energy levels are populated in the molecule's spectrum. However, it does NOT forbid the existence of the vibrational ground state that tunnels.
    """
    print(textwrap.dedent(explanation))

    # --- Final Conclusion ---
    print("\n--- Step 4: Final Conclusion ---")
    explanation = """
    Since the potential energy barrier for nitrogen inversion is still present in the exotic ammonia molecule, the fundamental condition for quantum tunneling is met. The Nitrogen atom's wavefunction can still penetrate this barrier.
    
    Therefore, the ammonia molecule with exotic, spin-0 hydrogens would still exhibit tunneling.
    """
    print(textwrap.dedent(explanation))


# Execute the explanation function
explain_ammonia_tunneling()