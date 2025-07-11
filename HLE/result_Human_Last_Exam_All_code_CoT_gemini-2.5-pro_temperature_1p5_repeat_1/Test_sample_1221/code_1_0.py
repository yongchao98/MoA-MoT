def calculate_residuary_estate():
    """
    Calculates the value of Betty's residuary estate based on the provided scenario.
    The calculation follows these steps:
    1. Sum up all assets belonging to Betty's estate at death.
    2. Sum up all valid legacies to be paid out from the will.
    3. Subtract total legacies from total assets to find the residuary estate.
    """

    print("--- Step 1: Calculating the Total Value of Betty's Estate Assets ---")

    # Flat A is owned as a joint tenancy. Under the commorientes rule in Hong Kong law
    # (Conveyancing and Property Ordinance, s.11), the younger person is presumed
    # to have survived the elder in an uncertain death scenario.
    # Betty (54) is younger than Alex (57), so she is presumed to have survived him.
    # The flat automatically passes to her by survivorship and becomes fully part of her estate.
    flat_a = 2000000
    print(f"Value of Flat A at Orchid Apartment Complex: {flat_a}")

    # Flat B is owned solely by Betty and is part of her estate.
    flat_b = 4000000
    print(f"Value of Flat B at Peach Apartment Complex: {flat_b}")

    # Other assets
    bank_accounts = 50000
    print(f"Value of cash in bank accounts: {bank_accounts}")
    nyse_shares = 30000
    print(f"Value of NYSE shares: {nyse_shares}")
    personal_items = 20000
    print(f"Value of personal items: {personal_items}")
    jewellery = 500000
    print(f"Value of jewellery collection: {jewellery}")

    # Calculate total assets
    total_assets = flat_a + flat_b + bank_accounts + nyse_shares + personal_items + jewellery
    print("\nCalculating Total Assets...")
    print(f"Total Assets = {flat_a} + {flat_b} + {bank_accounts} + {nyse_shares} + {personal_items} + {jewellery}")
    print(f"Total Assets = {total_assets}\n")

    print("--- Step 2: Calculating Total Valid Legacies to be Paid Out ---")

    # Clause 4: Gift of Flat C fails. The legal principle is 'ademption', as Betty no longer owned the specific item at her death.
    legacy_cindy_flat_c = 0
    print(f"Legacy of Flat C to Cindy: {legacy_cindy_flat_c} (Gift fails due to ademption as the flat was sold before death)")

    # Clause 5a: Gift to friends fails because the schedule listing the beneficiaries was never created. This fails for 'uncertainty of objects'.
    legacy_friends = 0
    print(f"Legacy to friends: {legacy_friends} (Gift fails for uncertainty as the schedule was not created)")

    # Clause 5b: This is a valid pecuniary legacy.
    legacy_wills_lawyers = 230000
    print(f"Legacy to Wills Lawyers & Co: {legacy_wills_lawyers} (This is a valid gift)")
    
    # Clause 6: This gift fails as the clause was deleted from the will.
    legacy_rspca = 0
    print(f"Legacy to RSPCA: {legacy_rspca} (Gift fails as the clause was deleted)")

    # Calculate total valid legacies
    total_legacies = legacy_cindy_flat_c + legacy_friends + legacy_wills_lawyers + legacy_rspca
    print("\nCalculating Total Valid Legacies...")
    print(f"Total Valid Legacies = {legacy_cindy_flat_c} + {legacy_friends} + {legacy_wills_lawyers} + {legacy_rspca}")
    print(f"Total Valid Legacies = {total_legacies}\n")
    
    print("--- Step 3: Calculating the Residuary Estate ---")
    print("The residuary estate is the amount remaining after all valid legacies, debts, and expenses are paid.")
    print("(Note: Debts are not quantified in the problem and are therefore excluded from this calculation.)")

    # Calculate the residuary estate
    residuary_estate = total_assets - total_legacies
    print("\nFinal Calculation:")
    print(f"Residuary Estate = Total Assets - Total Valid Legacies")
    print(f"Residuary Estate = {total_assets} - {total_legacies}")
    print(f"Final value of the residuary estate is: {residuary_estate}")

    print(f"\n<<<{residuary_estate}>>>")

# Run the calculation
calculate_residuary_estate()