import sys

def prove_instability():
    """
    This function uses a scaling argument (Derrick's Theorem) to prove that
    localized solitons are unstable in a 3D system with only Heisenberg
    exchange and DMI.
    """
    print("--- Proof of Soliton Instability via Scaling Argument ---")
    print("\nThe Hamiltonian density is given by: H = A*(grad m)^2 + D*m.(curl m)")
    print("The total energy is the volume integral of H: E = E_A + E_D\n")

    print("Step 1: Assume a localized soliton solution m_0(r) exists.")
    print("Let its exchange energy be E_A0 and its DMI energy be E_D0.")
    print("For any non-uniform soliton, the exchange energy must be positive.")
    print("Physical Condition: E_A0 = integral[A*(grad m_0)^2] dV > 0\n")

    print("Step 2: Introduce a scaled version of the solution: m_lambda(r) = m_0(r/lambda).")
    print("Here, 'lambda' is a positive scaling factor for the soliton's size.")
    print("We will now find the total energy E(lambda) of this scaled solution.\n")

    print("Step 3: Determine how the energy terms scale with lambda.")
    print("The gradient operator (grad) scales as 1/lambda.")
    print("The 3D volume element (dV) scales as lambda^3.\n")

    print("Scaling of Exchange Energy (E_A):")
    print("E_A(lambda) is proportional to (1/lambda)^2 for the gradient term and lambda^3 for the volume.")
    print("E_A(lambda) = E_A0 * (1/lambda)^2 * lambda^3")
    print("Resulting scaling: E_A(lambda) = lambda * E_A0\n")

    print("Scaling of DMI Energy (E_D):")
    print("E_D(lambda) is proportional to (1/lambda) for the curl term and lambda^3 for the volume.")
    print("E_D(lambda) = E_D0 * (1/lambda) * lambda^3")
    print("Resulting scaling: E_D(lambda) = lambda^2 * E_D0\n")

    print("Step 4: Formulate the total energy E(lambda) and find the equilibrium condition.")
    print("The total energy of the scaled soliton is:")
    print("E(lambda) = lambda * E_A0 + lambda^2 * E_D0\n")

    print("For the original soliton (lambda=1) to be an equilibrium state, the energy must be stationary.")
    print("This means the first derivative of E(lambda) must be zero at lambda=1.")
    print("dE/d(lambda) = E_A0 + 2 * lambda * E_D0")
    print("Setting dE/d(lambda) = 0 at lambda=1 gives the equilibrium condition (Virial Theorem):")
    print("Equilibrium Condition: E_A0 + 2 * E_D0 = 0\n")

    print("Step 5: Check the stability of this equilibrium.")
    print("For the soliton to be stable, the equilibrium must be a local minimum.")
    print("This means the second derivative of E(lambda) must be positive at lambda=1.")
    print("d^2E/d(lambda)^2 = 2 * E_D0")
    print("Stability Condition: 2 * E_D0 > 0  =>  E_D0 > 0\n")

    print("Step 6: Synthesize the results and identify the contradiction.")
    print("We have three conditions that must be met for a stable soliton:")
    print(f"1. Physical Condition: E_A0 > 0")
    print(f"2. Equilibrium Condition: E_A0 = -2 * E_D0")
    print(f"3. Stability Condition: E_D0 > 0\n")

    print("Let's check for consistency:")
    print("From the Stability Condition (3), E_D0 must be positive.")
    print("If E_D0 is positive, we can substitute it into the Equilibrium Condition (2):")
    print("E_A0 = -2 * (a positive number)")
    print("This implies that E_A0 must be negative.\n")

    print("--- CONCLUSION ---")
    print("The requirement that E_A0 is negative contradicts the Physical Condition (1) that E_A0 must be positive.")
    print("It is impossible to satisfy all conditions simultaneously.")
    print("Therefore, no stable localized soliton can exist in this system.")

if __name__ == '__main__':
    prove_instability()