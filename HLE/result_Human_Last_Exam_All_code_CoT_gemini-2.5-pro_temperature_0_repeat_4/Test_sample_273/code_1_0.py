import numpy as np

def analyze_soliton_stability():
    """
    Analyzes the stability of a 3D soliton based on Derrick's theorem.
    
    The analysis uses a scaling argument for the energy functional:
    E(lambda) = E_exchange(lambda) + E_DMI(lambda)
    
    In 3D, the terms scale with a characteristic size 'lambda' as:
    E_exchange(lambda) = E_ex_1 * lambda
    E_DMI(lambda) = E_dmi_1 * lambda^2
    """
    
    print("--- Soliton Stability Analysis in 3D ---")
    print("Hamiltonian: A(grad m)^2 + D(m . curl m)")
    print("\nStep 1: Define energy scaling with soliton size 'lambda'.")
    print("E_exchange(lambda) = E_ex_1 * lambda")
    print("E_DMI(lambda) = E_dmi_1 * lambda^2")
    print("Total Energy E(lambda) = (E_ex_1 * lambda) + (E_dmi_1 * lambda^2)")

    # For a soliton to form, the exchange energy must be positive (cost of gradients)
    # and the DMI energy must be negative (driving force for twisting).
    # Let's assume a hypothetical stationary solution exists at lambda=1.
    # The condition dE/d(lambda) = 0 at lambda=1 gives E_ex_1 + 2*E_dmi_1 = 0.
    # We can choose arbitrary values that satisfy this, e.g., E_ex_1 = 2, E_dmi_1 = -1.
    
    E_ex_1 = 2.0
    E_dmi_1 = -1.0
    
    print("\nStep 2: Check conditions for a stable minimum at lambda = 1.")
    print(f"Assuming a stationary point at lambda=1, we set E_ex_1 = {E_ex_1} and E_dmi_1 = {E_dmi_1}.")
    
    # First derivative of E(lambda) is E_ex_1 + 2 * E_dmi_1 * lambda
    dE_dlambda_at_1 = E_ex_1 + 2 * E_dmi_1 * 1.0
    print(f"First derivative at lambda=1: dE/d(lambda) = {E_ex_1} + 2*({E_dmi_1})*1 = {dE_dlambda_at_1}")
    if dE_dlambda_at_1 == 0:
        print(" -> Condition for a stationary point is met.")
    else:
        print(" -> Condition for a stationary point is NOT met.")

    # Second derivative of E(lambda) is 2 * E_dmi_1
    d2E_dlambda2_at_1 = 2 * E_dmi_1
    print(f"Second derivative at lambda=1: d^2E/d(lambda)^2 = 2*({E_dmi_1}) = {d2E_dlambda2_at_1}")
    if d2E_dlambda2_at_1 > 0:
        print(" -> The point is a local MINIMUM (stable).")
    else:
        print(" -> The point is a local MAXIMUM (unstable).")

    print("\nStep 3: Demonstrate the instability by calculating energy at different sizes.")
    
    def total_energy(l, e_ex, e_dmi):
        return e_ex * l + e_dmi * l**2

    # Energy at the stationary point (lambda=1)
    E_at_1 = total_energy(1.0, E_ex_1, E_dmi_1)
    print(f"Energy at lambda = 1.0 (hypothetical soliton): E = {E_at_1}")

    # Energy when shrinking (lambda -> 0)
    lambda_small = 0.1
    E_small = total_energy(lambda_small, E_ex_1, E_dmi_1)
    print(f"Energy at lambda = {lambda_small} (shrunken soliton): E = {E_small:.4f}. This is lower, indicating collapse instability.")

    # Energy when expanding (lambda -> infinity)
    lambda_large = 5.0
    E_large = total_energy(lambda_large, E_ex_1, E_dmi_1)
    print(f"Energy at lambda = {lambda_large} (expanded soliton): E = {E_large:.4f}. This is lower, indicating expansion instability.")

    print("\n--- Conclusion ---")
    print("The analysis shows that any localized soliton is at an energy maximum, not a minimum.")
    print("It is unstable and can lower its energy by either shrinking or expanding.")
    print("Therefore, it is not possible to stabilize such a soliton with only these two terms.")

if __name__ == '__main__':
    analyze_soliton_stability()