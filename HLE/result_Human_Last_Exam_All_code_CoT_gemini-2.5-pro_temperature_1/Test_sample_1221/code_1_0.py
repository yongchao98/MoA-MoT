def calculate_residuary_estate():
    """
    Calculates Betty's residuary estate based on the provided details.
    """
    # Step 1: Define and sum the assets in Betty's estate.
    # Flat A passes to Betty by survivorship (commorientes rule) as she was younger than Alex.
    flat_a = 2000000
    # Flat B was solely owned by Betty.
    flat_b = 4000000
    # Other assets.
    cash = 50000
    shares = 30000
    personal_items = 20000
    jewellery = 500000
    
    total_assets = flat_a + flat_b + cash + shares + personal_items + jewellery
    
    # Step 2: Identify valid legacies to be paid out from the estate.
    # Legacy for Flat C failed due to ademption (asset was sold).
    # Legacy to friends failed for uncertainty (schedule was never created).
    # Legacy to RSPCA was revoked (clause was deleted).
    
    # The only valid legacy is to Wills Lawyers & Co.
    legacy_wills_lawyers = 230000
    
    total_valid_legacies = legacy_wills_lawyers
    
    # Step 3: Calculate the residuary estate.
    residuary_estate = total_assets - total_valid_legacies
    
    # Print the breakdown of the calculation.
    print("Calculating Betty's Residuary Estate:")
    print("-" * 40)
    print("Assets included in the estate:")
    print(f"Flat A (by survivorship): {flat_a:,} HKD")
    print(f"Flat B: {flat_b:,} HKD")
    print(f"Cash: {cash:,} HKD")
    print(f"Shares: {shares:,} HKD")
    print(f"Personal Items: {personal_items:,} HKD")
    print(f"Jewellery: {jewellery:,} HKD")
    print(f"Total Estate Value: {total_assets:,} HKD")
    print("-" * 40)
    print("Deductions from the estate:")
    print(f"Valid Legacy to Wills Lawyers & Co: {total_valid_legacies:,} HKD")
    print("-" * 40)
    print("Final Equation for Residuary Estate:")
    print(f"{flat_a:,} + {flat_b:,} + {cash:,} + {shares:,} + {personal_items:,} + {jewellery:,} - {total_valid_legacies:,} = {residuary_estate:,}")
    print(f"\nThe total value of Betty's residuary estate is {residuary_estate:,} HKD.")

calculate_residuary_estate()
<<<6370000>>>