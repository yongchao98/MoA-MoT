import sys

def solve_lindhard_function():
    """
    This function explains and calculates the value of the Lindhard polarization
    function Π(q, ω) at q=0 and ω=0.
    """

    # The momentum transfer 'q' and frequency 'ω' are both zero.
    q_value = 0
    omega_value = 0

    # Based on physical reasoning, the result is 0.
    result = 0

    # Print the explanation of the result.
    print("In the context of linear response theory, the Lindhard polarization function, Π(q, ω), relates the induced charge density to an external potential.")
    print("\nWe are asked to find its value for zero momentum transfer (q=0) and zero frequency (ω=0).")
    print("\nPhysical Analysis:")
    print("1. A potential with q=0 and ω=0 is a uniform, static potential that adds a constant energy, U₀, to every electron.")
    print("2. To preserve the total number of electrons, the system's chemical potential, μ, must also shift by U₀.")
    print("3. The occupation of any state 'k' depends on the energy relative to the chemical potential, which is (ε_k + U₀) - (μ + U₀) = ε_k - μ.")
    print("4. Since this relative energy is unchanged, the occupation of every state is unchanged.")
    print("5. Therefore, the induced charge density is zero.")
    print("\nDeriving the value:")
    print("From the definition n_ind = Π(q, ω) * U_ext, we substitute our findings:")
    print(f"0 = Π(q={q_value}, ω={omega_value}) * U₀")
    print("For any non-zero potential U₀, this equation can only be true if Π(0, 0) is zero.")

    # Print the final result in an equation format.
    print("\nFinal Equation:")
    # sys.stdout.flush() # Ensure prints are out before the final result
    print(f"Π(q={q_value}, ω={omega_value}) = {result}")

solve_lindhard_function()