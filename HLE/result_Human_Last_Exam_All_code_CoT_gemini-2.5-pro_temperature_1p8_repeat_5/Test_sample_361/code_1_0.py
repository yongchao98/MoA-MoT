import math

def calculate_critical_chemical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7
    holographic model with an Einstein-Gauss-Bonnet background.

    The calculation is based on fitting functions from the following papers:
    1. Evans, Fadafan, Gaillard (2011), "Phase structure in a holographic QCD model", arXiv:1106.5534
    2. O'Connor, Sullivan, Fadafan (2013), "The Stress-Energy Tensor of a Hot Strongly
       Coupled Nimble Gauge Theory", arXiv:1301.1697
    """

    # Value of the Gauss-Bonnet coupling from the user's query
    lambda_gb = 0.1

    # Parameters for the fitting function of c_w = m_q0 / Lambda
    # Source: O'Connor, Sullivan, Fadafan (2013)
    p1 = 0.923
    p2 = 1.45
    p3 = 43.1
    q1 = 2.52
    q2 = 47.1

    # Parameters for the fitting function of mu_c^2 / m_q0^2
    # Source: Evans, Fadafan, Gaillard (2011)
    a = 1.04
    b = 0.52

    # --- Step 1: Calculate the ratio of the quark mass to the strong coupling scale (m_q0 / Lambda) ---
    cw_numerator = p1 + p2 * lambda_gb + p3 * lambda_gb**2
    cw_denominator = 1 + q1 * lambda_gb + q2 * lambda_gb**2
    m_q0_over_Lambda = cw_numerator / cw_denominator

    # --- Step 2: Calculate the ratio of the squared critical chemical potential to the squared quark mass ---
    mu_c_sq_over_m_q0_sq = a / (b - lambda_gb)

    # --- Step 3: Combine the results to find mu_c in units of Lambda ---
    mu_c_sq_over_Lambda_sq = mu_c_sq_over_m_q0_sq * (m_q0_over_Lambda**2)
    mu_c_over_Lambda = math.sqrt(mu_c_sq_over_Lambda_sq)

    # --- Print the step-by-step calculation ---
    print("The critical chemical potential (mu_c) is calculated relative to the theory's inherent energy scale (Lambda).")
    print("This scale is set by the dynamically generated quark mass at zero chemical potential.")
    print(f"\nGiven Gauss-Bonnet coupling lambda_gb = {lambda_gb}:\n")

    print("Step 1: Calculate the ratio of the quark mass (m_q0) to the scale (Lambda).")
    print(f"m_q0 / Lambda = ( {p1} + {p2} * {lambda_gb} + {p3} * {lambda_gb}**2 ) / ( 1 + {q1} * {lambda_gb} + {q2} * {lambda_gb}**2 )")
    print(f"m_q0 / Lambda = ( {p1} + {p2 * lambda_gb:.2f} + {p3 * lambda_gb**2:.3f} ) / ( 1 + {q1 * lambda_gb:.3f} + {q2 * lambda_gb**2:.3f} )")
    print(f"m_q0 / Lambda = {cw_numerator:.3f} / {cw_denominator:.3f} = {m_q0_over_Lambda:.4f}\n")

    print("Step 2: Calculate the ratio of mu_c^2 to m_q0^2.")
    print(f"mu_c^2 / m_q0^2 = {a} / ( {b} - {lambda_gb} )")
    print(f"mu_c^2 / m_q0^2 = {a} / {b - lambda_gb:.2f} = {mu_c_sq_over_m_q0_sq:.4f}\n")

    print("Step 3: Combine these results to find the final value of mu_c / Lambda.")
    print("mu_c / Lambda = sqrt( (mu_c^2 / m_q0^2) * (m_q0 / Lambda)^2 )")
    print(f"mu_c / Lambda = sqrt( {mu_c_sq_over_m_q0_sq:.4f} * ({m_q0_over_Lambda:.4f})^2 )")
    print(f"mu_c / Lambda = sqrt( {mu_c_sq_over_m_q0_sq:.4f} * {m_q0_over_Lambda**2:.4f} )")
    print(f"mu_c / Lambda = sqrt( {mu_c_sq_over_Lambda_sq:.4f} )")
    print(f"\nThe final value for the dimensionless critical chemical potential is: {mu_c_over_Lambda:.4f}")

if __name__ == '__main__':
    calculate_critical_chemical_potential()