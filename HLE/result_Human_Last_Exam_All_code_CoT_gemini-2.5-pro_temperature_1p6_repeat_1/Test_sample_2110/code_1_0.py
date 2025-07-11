import math

def calculate_cross_section(G_F, m_e, m_mu, s):
    """
    Calculates the total cross-section for e- + anti-nu_e -> mu- + anti-nu_mu.

    Args:
        G_F (float): Fermi constant.
        m_e (float): Mass of the electron.
        m_mu (float): Mass of the muon.
        s (float): Center-of-mass energy squared.
    """

    # Check if the process is kinematically allowed
    if s < m_mu**2:
        print("Center-of-mass energy s is below the muon production threshold.")
        return

    # Calculate energies and momenta in the CM frame
    # Incoming particles (electron, electron antineutrino)
    E_e = (s + m_e**2) / (2 * math.sqrt(s))
    p_in_mag = (s - m_e**2) / (2 * math.sqrt(s))
    
    # Outgoing particles (muon, muon antineutrino)
    E_mu = (s + m_mu**2) / (2 * math.sqrt(s))
    k_out_mag = (s - m_mu**2) / (2 * math.sqrt(s))

    # Calculate integration limits for Mandelstam variable t
    # t_max corresponds to cos(theta) = 1
    t_max_term = 2 * (E_e * E_mu - p_in_mag * k_out_mag)
    t_max = m_e**2 + m_mu**2 - t_max_term

    # t_min corresponds to cos(theta) = -1
    t_min_term = 2 * (E_e * E_mu + p_in_mag * k_out_mag)
    t_min = m_e**2 + m_mu**2 - t_min_term

    # The formula for the total cross-section is obtained by integrating d(sigma)/dt.
    # d(sigma)/dt = (1/(16*pi*(s-m_e^2)^2)) * <|M|^2>
    # <|M|^2> = -2*G_F^2 * t * (m_e^2 + m_mu^2 - t)
    # sigma = Integral from t_min to t_max of d(sigma)/dt
    #
    # Performing the integral:
    # I = Integral [ (m_e^2 + m_mu^2)*t - t^2 ] dt
    # I = [ (m_e^2 + m_mu^2)*(t^2/2) - t^3/3 ] from t_min to t_max

    C = m_e**2 + m_mu**2
    
    def F(t):
        return C * (t**2 / 2) - t**3 / 3
    
    integral_val = F(t_max) - F(t_min)
    
    sigma_factor = (-G_F**2) / (8 * math.pi * (s - m_e**2)**2)
    
    sigma = sigma_factor * integral_val

    # Display the parameters and the result
    print(f"Given parameters:")
    print(f"G_F = {G_F}")
    print(f"m_e = {m_e}")
    print(f"m_μ = {m_mu}")
    print(f"s = {s}")
    print("\nCalculation steps:")
    print(f"Center-of-mass energies and momenta:")
    print(f"E_e = {E_e:.4f}, E_μ = {E_mu:.4f}")
    print(f"|p_in| = {p_in_mag:.4f}, |k_out| = {k_out_mag:.4f}")
    print(f"Integration limits for t:")
    print(f"t_min = {t_min:.4f}")
    print(f"t_max = {t_max:.4f}")
    
    print("\nThe cross-section formula is:")
    print(f"σ = [-({G_F})² / (8π * ({s} - {m_e}²)²)] * [({m_e}²+{m_mu}²)t²/2 - t³/3] from t={t_min:.2f} to {t_max:.2f}")

    print(f"\nFinal result:")
    print(f"Total cross-section σ = {sigma:.6f}")
    
    return sigma

# Given values
G_F_val = 1
m_e_val = 1
m_mu_val = 1
s_val = 2

final_sigma = calculate_cross_section(G_F_val, m_e_val, m_mu_val, s_val)
print(f"\n<<<{final_sigma:.6f}>>>")