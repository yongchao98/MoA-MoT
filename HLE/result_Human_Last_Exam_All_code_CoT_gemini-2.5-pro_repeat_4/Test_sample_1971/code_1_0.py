import sympy as sp

def solve_sugra_parameters():
    """
    This function follows the plan to derive the values of alpha^2 and beta.
    It uses strings to explain the steps and print the final results,
    as requested by the user.
    """

    print("Step-by-step derivation of alpha^2 and beta:")
    print("-" * 50)

    # Step 1: Find the vacuum expectation value (VEV) of S
    print("1. Find the VEV of the auxiliary field S, denoted S_0.")
    print("The Lagrangian terms involving S are L(S) = -1/3 * e * S^2 + alpha * e * S.")
    print("The equation of motion for S is dL/dS = 0, which gives:")
    print("-2/3 * e * S_0 + alpha * e = 0")
    print("Solving for S_0, we get:")
    print("S_0 = (3/2) * alpha")
    S0_expr = "3/2 * alpha"
    print("-" * 50)

    # Step 2: Impose supersymmetry of the vacuum (delta_psi_mu = 0)
    print("2. Impose that the vacuum is supersymmetric, i.e., delta_psi_mu = 0.")
    print("The gravitino transformation is: delta_psi_mu = (1/kappa) * D_mu(epsilon) + 1/6 * gamma_mu * S * epsilon.")
    print("In the vacuum (psi_mu = 0, S = S_0), this becomes:")
    print("(1/kappa) * D_mu(epsilon) + 1/6 * gamma_mu * S_0 * epsilon = 0")
    print("This gives the Killing spinor equation: D_mu(epsilon) = - (kappa * S_0 / 6) * gamma_mu * epsilon.")
    print(f"Substituting S_0 = {S0_expr}:")
    print("D_mu(epsilon) = - (kappa * (3/2 * alpha) / 6) * gamma_mu * epsilon = - (kappa * alpha / 4) * gamma_mu * epsilon.")
    m_coeff_str = "-(kappa * alpha / 4)"
    print("-" * 50)

    # Step 3: Use the integrability condition to find alpha^2
    print("3. Use the integrability condition [D_mu, D_nu]epsilon = 1/2 * R_munu^ab * gamma_ab * epsilon.")
    print("The left-hand side gives: [D_mu, D_nu]epsilon = (-kappa*alpha/4)^2 * [gamma_nu, gamma_mu]epsilon")
    print("= (kappa^2 * alpha^2 / 16) * (-2 * gamma_munu) * epsilon = -(kappa^2 * alpha^2 / 8) * gamma_munu * epsilon.")
    print("Equating the two sides implies the spacetime has constant curvature.")
    print("For a space of constant curvature, the Ricci scalar R is related to the coefficient by R = -3 * kappa^2 * alpha^2.")
    print("Solving for alpha^2 gives:")
    alpha_squared_val_1 = -1
    alpha_squared_val_2 = 3
    print(f"alpha^2 = ({alpha_squared_val_1} * R) / ({alpha_squared_val_2} * kappa^2)")
    print("-" * 50)

    # Step 4: Determine beta from the gravitino mass term
    print("4. Determine beta by comparing the gravitino mass term in the Lagrangian.")
    print("From the Killing spinor equation, the gravitino mass m_3/2 is defined by D_mu(epsilon) = -(m_3/2 / 2) * gamma_mu * epsilon (by convention, taking the positive root).")
    print(f"Comparing with D_mu(epsilon) = {m_coeff_str} * gamma_mu * epsilon, we get m_3/2 = kappa * alpha / 2.")
    print("The canonical mass term in the Lagrangian is L_mass_canonical = -e/2 * m_3/2 * bar(psi)_mu * gamma^munu * psi_nu.")
    print("The mass term in our L_cos is L_mass_cos = alpha * e * kappa * beta * bar(psi)_mu * gamma^munu * psi_nu.")
    print("Equating the coefficients: alpha * kappa * beta = -m_3/2 / 2.")
    print("alpha * kappa * beta = -(kappa * alpha / 2) / 2 = -kappa * alpha / 4.")
    print("Solving for beta gives:")
    beta_val_1 = -1
    beta_val_2 = 4
    print(f"beta = {beta_val_1} / {beta_val_2}")
    print("-" * 50)

    print("\nFinal Results:")
    print("The value of alpha^2 in terms of the constant curvature R and kappa is:")
    print(f"alpha^2 = ({alpha_squared_val_1} * R) / ({alpha_squared_val_2} * kappa^2)")
    print(f"The numbers in the final equation for alpha^2 are: {alpha_squared_val_1}, {alpha_squared_val_2}")
    print("\nThe value of beta is:")
    print(f"beta = {beta_val_1} / {beta_val_2}")
    print(f"The numbers in the final equation for beta are: {beta_val_1}, {beta_val_2}")


if __name__ == "__main__":
    solve_sugra_parameters()
