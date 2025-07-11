import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_bose_einstein_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons
    in Bose-Einstein condensation and identifies the correct answer.
    """
    # The user is asking for the fundamental limit on the chemical potential (μ) for bosons in Bose-Einstein condensation (BEC).
    # This is a core concept in statistical mechanics, derived from the Bose-Einstein distribution.

    # Step 1: Explain the Bose-Einstein distribution.
    # The average number of bosons in a single-particle state 'i' with energy ϵ_i is given by the Bose-Einstein distribution.
    # We will print this equation and explain its components.
    print("In the grand canonical ensemble, the average number of bosons in a single-particle state 'i' with energy ϵ_i is given by the Bose-Einstein distribution:")
    print("n_i = 1 / (exp((ϵ_i - μ) / (k_B * T)) - 1)")
    print("where:")
    print("  μ is the chemical potential")
    print("  ϵ_i is the energy of the state 'i'")
    print("  k_B is the Boltzmann constant")
    print("  T is the temperature\n")

    # Step 2: Derive the constraint on μ from physical principles.
    # The occupation number n_i must be non-negative. This places a constraint on the denominator.
    print("For the occupation number n_i to be a positive physical value (n_i >= 0), the denominator of this expression must be positive:")
    print("exp((ϵ_i - μ) / (k_B * T)) - 1 > 0")
    print("This implies:")
    print("exp((ϵ_i - μ) / (k_B * T)) > 1\n")

    print("Taking the natural logarithm of both sides gives:")
    print("(ϵ_i - μ) / (k_B * T) > 0")
    print("Since both k_B and T are positive, this simplifies to the fundamental condition on the chemical potential:")
    print("μ < ϵ_i\n")

    # Step 3: Apply the constraint to the ground state.
    # This condition must hold for ALL energy states, including the one with the lowest energy (the ground state, ϵ_0).
    print("This inequality must hold for all possible energy states. Therefore, the chemical potential must be less than the lowest possible energy, the ground state energy ϵ_0:")
    print("μ < ϵ_0\n")

    # Step 4: Explain what happens at condensation.
    # For a macroscopic number of particles to occupy the ground state (n_0 -> ∞), the denominator for n_0 must approach zero.
    # This means μ must get very, very close to ϵ_0.
    print("For Bose-Einstein condensation to occur, a macroscopic number of particles must occupy the ground state. Let's look at the equation for n_0:")
    print("n_0 = 1 / (exp((ϵ_0 - μ) / (k_B * T)) - 1)")
    print("For n_0 to become very large, its denominator must approach zero. This happens as the chemical potential μ approaches the ground state energy ϵ_0 from below.")
    print("In the condensed phase (for temperatures T below the critical temperature T_c), the chemical potential is effectively pinned at this limit: μ = ϵ_0.\n")

    # Step 5: Connect to the answer choices.
    # The chemical potential of a non-interacting Bose gas at absolute zero (T=0) is exactly the ground state energy, ϵ_0.
    # Therefore, the condition μ = ϵ_0 is equivalent to saying μ equals the chemical potential of a non-interacting Bose gas at T=0.
    print("The chemical potential of a non-interacting Bose gas at absolute zero temperature (T=0) is precisely the ground state energy, μ(T=0) = ϵ_0.")
    print("Therefore, the limit that defines the condensed phase is that the chemical potential μ becomes equal to the chemical potential of a non-interacting Bose gas at zero temperature.")
    print("This corresponds to choice C.")


solve_bose_einstein_limit()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the explanation and the final answer in the required format.
print(output)
print("<<<C>>>")