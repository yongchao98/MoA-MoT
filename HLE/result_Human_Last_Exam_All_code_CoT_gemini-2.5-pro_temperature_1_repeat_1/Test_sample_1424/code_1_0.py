import cmath

def calculate_bare_greens_function(omega, epsilon_k, eta):
    """
    Calculates the bare Green's function G_0(k, omega).

    Args:
        omega (float): The energy variable.
        epsilon_k (float): The single-particle energy eigenvalue.
        eta (float): A small positive infinitesimal for causality.
    """
    # The formula is G_0 = 1 / (omega - epsilon_k + i*eta)
    denominator = (omega - epsilon_k) + 1j * eta
    g0 = 1 / denominator

    print("The functional dependence of the bare Green's function G_0 on the single-particle energy")
    print("eigenvalue epsilon_k is given by:")
    print("G_0 = 1 / (omega - epsilon_k + i*eta)\n")

    print(f"Using the following example values:")
    print(f"omega = {omega}")
    print(f"epsilon_k = {epsilon_k}")
    print(f"eta = {eta}\n")

    # Print the equation with the specific numbers plugged in
    print("The calculation is:")
    print(f"G_0 = 1 / (({omega}) - ({epsilon_k}) + i*({eta}))")
    print(f"G_0 = 1 / ({denominator:.4f})")
    print(f"G_0 = {g0:.4f}")

# --- Parameters ---
# Energy variable (e.g., in eV)
omega = 2.5
# Single-particle energy eigenvalue (e.g., in eV)
epsilon_k = 2.2
# Infinitesimal shift (a small positive number)
eta = 0.1

# --- Calculation and Output ---
calculate_bare_greens_function(omega, epsilon_k, eta)
