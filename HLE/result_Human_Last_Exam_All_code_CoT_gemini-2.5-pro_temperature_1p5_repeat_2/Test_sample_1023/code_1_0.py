def solve_diversification_rates():
    """
    Calculates the multiplicative factor by which the morphospecies extinction rate
    is greater than the evolutionary species extinction rate based on the problem's context.
    """

    # Let mu_e be the normalized rate of true extinction for an evolutionary species. We can set it to 1.
    mu_e_rate = 1

    # Let mu_m be the extinction rate for a morphospecies.
    # mu_m = (rate of true extinction) + (rate of anagenetic pseudo-extinction) + (rate of bifurcating pseudo-extinction)
    # mu_m = mu_e + sigma + 0.5 * lambda_e

    # The total rate of pseudo-extinction is (sigma + 0.5 * lambda_e).
    # The key insight is to assume that the rate of pseudo-extinction (an artifact of classification)
    # is equal to the rate of true extinction.
    # Assumption: pseudo_extinction_rate = mu_e_rate
    pseudo_extinction_rate = mu_e_rate
    
    # Now, we calculate the total morphospecies extinction rate based on this assumption.
    # mu_m = mu_e + (pseudo_extinction_rate)
    mu_m_rate = mu_e_rate + pseudo_extinction_rate

    # The multiplicative factor is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.
    factor = mu_m_rate / mu_e_rate

    print("Step 1: The extinction rate for an evolutionary species is the true extinction rate, mu_e.")
    print(f"Let's represent the normalized rate of mu_e as: {mu_e_rate}")
    print("\nStep 2: The extinction rate for a morphospecies (mu_m) is the sum of true extinction (mu_e) and pseudo-extinction rates.")
    print("Pseudo-extinction comes from anagenesis (rate sigma) and bifurcation (rate 0.5*lambda_e).")
    print("Equation: mu_m = mu_e + (sigma + 0.5*lambda_e)")

    print("\nStep 3: Assume the total pseudo-extinction rate equals the true extinction rate.")
    print(f"Assumption: (sigma + 0.5*lambda_e) = mu_e = {pseudo_extinction_rate}")
    
    print("\nStep 4: Substitute the assumption into the equation for mu_m.")
    print(f"mu_m = mu_e + mu_e = {mu_e_rate} + {pseudo_extinction_rate} = {int(mu_m_rate)}")

    print("\nStep 5: Calculate the final factor (mu_m / mu_e).")
    print(f"Factor = {int(mu_m_rate)} / {mu_e_rate} = {int(factor)}")


solve_diversification_rates()
<<<2>>>