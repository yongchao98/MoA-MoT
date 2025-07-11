import math

def calculate_nsvz_beta_function(g, Nc, Nf, gamma):
    """
    Calculates the NSVZ beta function for an SU(Nc) gauge theory with Nf matter flavors.

    Args:
        g (float): The gauge coupling constant.
        Nc (int): The number of colors (e.g., 3 for SU(3)).
        Nf (int): The number of matter flavors in the fundamental representation.
        gamma (float): The anomalous dimension of the matter superfields.
    """
    pi = math.pi
    
    # Dynkin index for the adjoint representation of SU(Nc)
    T_G = Nc
    
    # Dynkin index for the fundamental representation of SU(Nc)
    T_R = 0.5
    
    # The NSVZ formula is:
    # beta(g) = - (g^3 / (16*pi^2)) * (Numerator / Denominator)
    # Numerator = 3*T(G) - sum_i[T(R_i) * (1 - gamma_i)]
    # Denominator = 1 - T(G) * g^2 / (8*pi^2)
    
    print("--- Calculating NSVZ Beta Function ---")
    print(f"Parameters: Gauge Coupling g = {g}, Nc = {Nc}, Nf = {Nf}, Anomalous Dimension gamma = {gamma}")
    
    # Calculate the numerator part of the formula
    # We assume all Nf flavors have the same anomalous dimension gamma
    numerator = 3 * T_G - Nf * T_R * (1 - gamma)
    
    # Calculate the denominator part of the formula
    denominator = 1 - (T_G * g**2) / (8 * pi**2)
    
    # Calculate the overall prefactor
    prefactor = - (g**3) / (16 * pi**2)
    
    # Calculate the final beta function value
    beta_g = prefactor * (numerator / denominator)
    
    print("\n--- Equation Components ---")
    print(f"Prefactor -(g^3 / 16*pi^2) = {prefactor:.6f}")
    print(f"Numerator (3*T(G) - Nf*T(R)*(1-gamma)) = {numerator:.4f}")
    print(f"Denominator (1 - T(G)*g^2 / (8*pi^2)) = {denominator:.6f}")
    print("\n--- Final Result ---")
    print(f"The value of the beta function, beta(g), is: {beta_g:.6f}")

# --- Example Usage ---
# Let's consider a theory similar to Supersymmetric QCD
if __name__ == '__main__':
    gauge_coupling = 1.2  # Example value for the gauge coupling
    num_colors = 3        # SU(3)
    num_flavors = 5       # 5 flavors of quarks/squarks
    anomalous_dim = -0.1  # Example value for the anomalous dimension
    
    calculate_nsvz_beta_function(gauge_coupling, num_colors, num_flavors, anomalous_dim)
