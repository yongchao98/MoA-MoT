def model_aldh_change():
    """
    This script models the change in ALDH levels in response to two different
    electrophilic compounds that activate the Keap1-Nrf2 pathway.
    """
    # --- Initial Parameters ---
    # Baseline ALDH level in arbitrary units
    base_aldh = 100.0
    # Concentration of treatment in uM
    concentration = 50.0

    # --- Potency Factors (Conceptual) ---
    # Based on biological literature, 4-OI is a more potent Nrf2 activator
    # than HNE derivatives. These factors represent relative potency per uM.
    potency_hne_alkyne = 0.025
    potency_4_oi = 0.040

    # --- Calculation for (2E)-4-Hydroxy-2-nonen-8-ynal ---
    change_hne_alkyne = base_aldh * potency_hne_alkyne * concentration
    final_aldh_hne_alkyne = base_aldh + change_hne_alkyne

    # --- Calculation for 4-OI ---
    change_4_oi = base_aldh * potency_4_oi * concentration
    final_aldh_4_oi = base_aldh + change_4_oi

    # --- Output Results ---
    print("Modeling ALDH change via the Keap1-Nrf2 pathway.")
    print("="*50)
    print("1. Treatment with (2E)-4-Hydroxy-2-nonen-8-ynal:")
    print(f"   - This electrophile activates Nrf2, causing an INCREASE in ALDH.")
    print(f"   - Modeled Equation: Final ALDH = Base ALDH + (Base ALDH * Potency * Concentration)")
    print(f"   - Calculation: {final_aldh_hne_alkyne:.1f} = {base_aldh:.1f} + ({base_aldh:.1f} * {potency_hne_alkyne} * {concentration:.1f})")
    print(f"   - Resulting increase: {change_hne_alkyne:.1f} units")
    print("-"*50)

    print("2. Treatment with 4-OI:")
    print(f"   - 4-OI is a MORE potent Nrf2 activator.")
    print(f"   - Modeled Equation: Final ALDH = Base ALDH + (Base ALDH * Potency * Concentration)")
    print(f"   - Calculation: {final_aldh_4_oi:.1f} = {base_aldh:.1f} + ({base_aldh:.1f} * {potency_4_oi} * {concentration:.1f})")
    print(f"   - Resulting increase: {change_4_oi:.1f} units")
    print("="*50)

    print("\nConclusion:")
    print(f"The change in ALDH is an INCREASE for both compounds.")
    print(f"The increase from 4-OI ({change_4_oi:.1f}) is MORE than from HNE-alkyne ({change_hne_alkyne:.1f}).")
    print(f"The key protein sensor in this pathway is Keap1.")
    print("\nFinal Answer components: increase, more, Keap1")

if __name__ == '__main__':
    model_aldh_change()
<<<B>>>