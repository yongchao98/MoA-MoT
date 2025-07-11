import math

def calculate_critical_chemical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model using an analytical approximation.
    
    The formula used is from S. Kanno, Phys. Rev. D 84, 086004 (2011),
    derived using a variational method. Please note that this is an
    approximation; numerical shooting methods yield a more accurate value of
    approximately 3.65 for the same parameters.
    """
    
    # --- Input Parameters ---
    # Gauss-Bonnet coupling
    lambda_gb = 0.1
    # AdS radius (set to 1)
    L = 1.0
    # Conformal dimension of the condensing operator (e.g., <q-bar q> has Delta=3)
    Delta = 3
    # Boundary spacetime dimension
    d = 4
    
    # --- Step 1: Calculate the effective AdS radius L_eff ---
    # The relation is L_eff^2 = (2 * lambda_gb * L^2) / (1 - sqrt(1 - 4 * lambda_gb))
    if 1 - 4 * lambda_gb < 0:
        print("Error: The Gauss-Bonnet coupling is too large. The AdS background is not well-defined.")
        return
        
    l_eff_sq = (2 * lambda_gb * L**2) / (1 - math.sqrt(1 - 4 * lambda_gb))
    l_eff = math.sqrt(l_eff_sq)
    
    # --- Step 2: Calculate the squared mass of the scalar field ---
    # The mass is related to the operator dimension by m^2 * L^2 = Delta * (Delta - d)
    m_sq = Delta * (Delta - d) / L**2
    
    # --- Step 3: Calculate the effective dimension Delta_plus ---
    # Delta_plus is given by 2 + sqrt(4 + m^2 * L_eff^2)
    m_sq_l_eff_sq = m_sq * l_eff_sq
    if 4 + m_sq_l_eff_sq < 0:
        print("Error: The effective mass squared violates the BF bound.")
        return
        
    delta_plus = 2 + math.sqrt(4 + m_sq_l_eff_sq)

    # --- Step 4: Calculate the critical chemical potential mu_c ---
    # The formula is mu_c = (Delta_+ * (Delta_+ - 1)) / (sqrt(2*Delta_*(Delta_+-1) - Delta_+) * L_eff)
    
    term1 = delta_plus * (delta_plus - 1)
    
    inner_sqrt_term = 2 * term1 - delta_plus
    if inner_sqrt_term < 0:
        print("Error: Cannot compute mu_c due to negative square root.")
        return
    
    denominator = math.sqrt(inner_sqrt_term) * l_eff
    
    mu_c = term1 / denominator
    
    # --- Step 5: Output the final equation and result ---
    print("The critical chemical potential mu_c is calculated using the equation:")
    print(f"mu_c = ({delta_plus:.3f} * ({delta_plus:.3f} - 1)) / (sqrt(2 * {delta_plus:.3f} * ({delta_plus:.3f} - 1) - {delta_plus:.3f}) * {l_eff:.3f})")
    print("\nPlugging in the numbers:")
    print(f"mu_c = ({term1:.3f}) / (sqrt({2*term1:.3f} - {delta_plus:.3f}) * {l_eff:.3f})")
    print(f"mu_c = {term1:.3f} / ({math.sqrt(inner_sqrt_term):.3f} * {l_eff:.3f})")
    print(f"mu_c = {term1:.3f} / {denominator:.3f}")
    print(f"\nThe calculated critical chemical potential is: {mu_c:.3f}")

if __name__ == '__main__':
    calculate_critical_chemical_potential()