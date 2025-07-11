import numpy as np

def calculate_bare_greens_function(omega, epsilon_k, mu, eta):
    """
    Calculates the bare Green's function G_0(k, omega).

    Args:
        omega (float): The probe particle energy.
        epsilon_k (float): The single-particle energy eigenvalue.
        mu (float): The chemical potential.
        eta (float): A small positive infinitesimal.

    Returns:
        complex: The value of the Green's function.
    """
    energy_diff = epsilon_k - mu
    sign = np.sign(energy_diff)
    
    # Handle the case where epsilon_k == mu, where the sign is 0
    if sign == 0:
        # At the Fermi level, the pole structure is more complex.
        # For this demonstration, we can treat it as the unoccupied case,
        # but in reality it's a mix of both. We'll set sign to 1.
        sign = 1

    denominator = omega - energy_diff + 1j * eta * sign
    g0 = 1 / denominator
    return g0, energy_diff, sign, denominator

def print_calculation(omega, epsilon_k, mu, eta):
    """Prints a step-by-step calculation of the Green's function."""
    
    g0, energy_diff, sign, denominator = calculate_bare_greens_function(omega, epsilon_k, mu, eta)

    state_type = "unoccupied" if energy_diff > 0 else "occupied"
    print(f"--- Calculating for an {state_type} state ---")
    print(f"The formula is: G_0(k, w) = 1 / (w - (e_k - mu) + i*eta*sign(e_k - mu))")
    print("\nSubstituting the given values:")
    print(f"w = {omega}")
    print(f"e_k = {epsilon_k}")
    print(f"mu = {mu}")
    print(f"eta = {eta}")
    
    print("\nStep-by-step calculation:")
    print(f"1. Energy relative to chemical potential: (e_k - mu) = ({epsilon_k} - {mu}) = {energy_diff:.2f}")
    print(f"2. Sign of (e_k - mu): sign({energy_diff:.2f}) = {int(sign)}")
    print(f"3. Denominator: ({omega} - {energy_diff:.2f}) + 1j * {eta} * {int(sign)}")
    print(f"   = {denominator.real:.4f} + {denominator.imag:.4f}j")
    print(f"4. Final Result: G_0 = 1 / ({denominator.real:.4f} + {denominator.imag:.4f}j)")
    print(f"   G_0 = {g0.real:.4f} + {g0.imag:.4f}j")
    print("-" * 35 + "\n")


if __name__ == "__main__":
    # Define parameters for the demonstration
    omega_val = 1.0  # Probe energy
    mu_val = 0.5   # Chemical potential
    eta_val = 0.1  # Small infinitesimal (made larger for clear demonstration)

    # --- Case 1: Unoccupied State (e_k > mu) ---
    epsilon_k_unoccupied = 1.2
    print_calculation(omega_val, epsilon_k_unoccupied, mu_val, eta_val)

    # --- Case 2: Occupied State (e_k < mu) ---
    epsilon_k_occupied = -0.2
    print_calculation(omega_val, epsilon_k_occupied, mu_val, eta_val)