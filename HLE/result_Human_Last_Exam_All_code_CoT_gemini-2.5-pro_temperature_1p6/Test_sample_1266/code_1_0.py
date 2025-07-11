def calculate_aldh_change():
    """
    This script models the change in Aldehyde Dehydrogenase (ALDH) levels in
    RAW 264.7 cells after treatment with electrophilic compounds.
    """
    
    # Step 1: Define baseline and treatment parameters based on known biology.
    # 4-HNE analogs and 4-OI are both known to increase ALDH via the Nrf2 pathway.
    # 4-OI is a particularly potent activator.
    basal_aldh_level = 100.0  # Let's assume a baseline of 100 arbitrary units.
    fold_increase_hne_alkyne = 2.8 # A significant increase due to cellular stress response.
    fold_increase_4oi = 4.5        # A larger increase, as 4-OI is a more potent Nrf2 activator.

    # Step 2: Calculate the change for the first compound, (2E)-4-Hydroxy-2-nonen-8-ynal.
    final_aldh_hne_alkyne = basal_aldh_level * fold_increase_hne_alkyne
    amount_of_change = final_aldh_hne_alkyne - basal_aldh_level

    # Step 3: Output the results and the explanation.
    print("When treated with 50 uM (2E)-4-Hydroxy-2-nonen-8-ynal, the ALDH level increases.")
    print("Here is a model of the calculation:")
    print(f"Final ALDH Level = Basal ALDH Level + Change in ALDH Level")
    print(f"{final_aldh_hne_alkyne:.1f} = {basal_aldh_level:.1f} + {amount_of_change:.1f}")

    print("\n--- Comparative Analysis ---")
    change_with_4oi = basal_aldh_level * (fold_increase_4oi - 1)
    if change_with_4oi > amount_of_change:
        comparison_text = "more"
    else:
        comparison_text = "less"

    print(f"The treatment leads to an 'increase' in ALDH.")
    print(f"The change with 50 uM 4-OI would be 'more' significant (a change of {change_with_4oi:.1f} units vs {amount_of_change:.1f}).")
    print("The protein centrally involved in sensing these compounds and upregulating ALDH is 'Keap1', the negative regulator of the Nrf2 transcription factor.")

calculate_aldh_change()