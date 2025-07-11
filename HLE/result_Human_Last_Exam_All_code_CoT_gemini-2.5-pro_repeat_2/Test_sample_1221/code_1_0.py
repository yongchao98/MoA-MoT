def calculate_residuary_estate():
    """
    Calculates the value of Betty's residuary estate based on the provided facts and Hong Kong law.
    """
    # Step 1: Determine the total value of assets in Betty's estate (in HKD).

    # Flat A: Owned as joint tenants. Betty (54) is presumed to have survived Alex (57).
    # The right of survivorship applies, so the full value of the flat goes into Betty's estate.
    flat_a = 2000000

    # Flat B: Owned solely by Betty.
    flat_b = 4000000

    # Other assets
    bank_accounts = 50000
    shares = 30000
    personal_items = 20000
    jewellery = 500000

    # Sum of all assets
    total_assets = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery

    # Step 2: Determine the value of valid legacies to be paid out.

    # Gift of Flat C (Clause 4): Fails due to ademption as the flat was sold. Value = 0.
    gift_cindy = 0

    # Gift to friends (Clause 5a): Fails for uncertainty of objects as the schedule was not created. Value = 0.
    gift_friends = 0

    # Gift to Wills Lawyers & Co (Clause 5b): This is a valid pecuniary legacy.
    gift_lawyers = 230000

    # Gift to RSPCA (Clause 6): The clause was deleted and is therefore revoked. Value = 0.
    gift_rspca = 0

    # Sum of all valid legacies
    total_legacies = gift_cindy + gift_friends + gift_lawyers + gift_rspca

    # Step 3: Calculate the residuary estate.
    residuary_estate = total_assets - total_legacies
    
    # Print the breakdown of the calculation as requested.
    print("Calculating Betty's Residuary Estate (in HKD)")
    print("---------------------------------------------")
    print("Assets:")
    print(f"  Flat A (by survivorship): {flat_a}")
    print(f"  Flat B:                   {flat_b}")
    print(f"  Bank Accounts:            {bank_accounts}")
    print(f"  Shares:                   {shares}")
    print(f"  Personal Items:           {personal_items}")
    print(f"  Jewellery:                {jewellery}")
    print(f"Total Assets = {flat_a} + {flat_b} + {bank_accounts} + {shares} + {personal_items} + {jewellery} = {total_assets}")
    print("\nLegacies Paid Out:")
    print(f"  Legacy to Wills Lawyers & Co: {gift_lawyers}")
    print(f"Total Legacies Paid = {total_legacies}")
    print("\nFinal Calculation:")
    print(f"Residuary Estate = Total Assets - Total Legacies Paid")
    print(f"Final Equation: {total_assets} - {total_legacies} = {residuary_estate}")
    print(f"\nThe total value of Betty's residuary estate is {residuary_estate} HKD.")


calculate_residuary_estate()
<<<6370000>>>