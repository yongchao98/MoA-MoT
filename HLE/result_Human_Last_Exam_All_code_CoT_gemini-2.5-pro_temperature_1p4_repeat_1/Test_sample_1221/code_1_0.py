def calculate_residuary_estate():
    """
    Calculates the value of Betty's residuary estate based on the provided details.
    """

    # --- Step 1 & 2: Identify and value assets forming the estate ---

    # Betty (age 54) is younger than Alex (age 57). Under the commorientes rule
    # (s.11 of the Conveyancing and Property Ordinance), Betty is presumed to have survived Alex.
    # Therefore, the jointly owned Flat A passes to Betty by survivorship and becomes part of her estate.
    flat_a_orchid_apartment = 2000000

    # Betty's solely owned assets
    flat_b_peach_apartment = 4000000
    bank_accounts = 50000
    nyse_shares = 30000
    personal_items = 20000
    jewellery_collection = 500000

    # --- Step 3: Sum the value of all assets in the estate ---
    total_assets = (
        flat_a_orchid_apartment +
        flat_b_peach_apartment +
        bank_accounts +
        nyse_shares +
        personal_items +
        jewellery_collection
    )

    # --- Step 4: Analyze gifts in the will ---

    # Clause 4: Gift of Flat C to Cindy Tran. This gift is adeemed (fails) because the property was sold.
    gift_to_cindy = 0

    # Clause 5a: Gift to friends. This gift is void for uncertainty as the schedule was not created.
    gift_to_friends = 0

    # Clause 5b: Gift to Wills Lawyers & Co. This is a valid pecuniary legacy.
    gift_to_lawyers = 230000

    # Clause 6: Gift to RSPCA. This clause was deleted from the will.
    gift_to_rspca = 0

    total_valid_legacies = gift_to_cindy + gift_to_friends + gift_to_lawyers + gift_to_rspca

    # --- Step 5: Calculate the residuary estate ---
    residuary_estate = total_assets - total_valid_legacies

    # --- Print the results ---
    print("Calculation of Betty's Residuary Estate:\n")

    print("Assets Forming the Estate:")
    print(f"- Flat A (passes to Betty's estate by survivorship): {flat_a_orchid_apartment:10,d} HKD")
    print(f"- Flat B (solely owned):                         {flat_b_peach_apartment:10,d} HKD")
    print(f"- Bank Accounts:                                  {bank_accounts:10,d} HKD")
    print(f"- NYSE Shares:                                    {nyse_shares:10,d} HKD")
    print(f"- Personal Items:                                 {personal_items:10,d} HKD")
    print(f"- Jewellery Collection:                           {jewellery_collection:10,d} HKD")
    print("-------------------------------------------------")
    print(f"Total Estate Value:                               {total_assets:10,d} HKD\n")

    print("Deductions for Legacies from Will:")
    print(f"- Legacy to Wills Lawyers & Co:                   {total_valid_legacies:10,d} HKD")
    print("  (Other legacies failed due to ademption, uncertainty, or deletion)\n")

    print("Final Residuary Estate Calculation:")
    print(f"Equation: {flat_a_orchid_apartment} + {flat_b_peach_apartment} + {bank_accounts} + {nyse_shares} + {personal_items} + {jewellery_collection} - {total_valid_legacies}")
    print(f"Result:   {total_assets:,} - {total_valid_legacies:,} = {residuary_estate:,} HKD")
    print("\n-------------------------------------------------")
    print(f"The total residuary estate is {residuary_estate:,} HKD.")


if __name__ == "__main__":
    calculate_residuary_estate()
    # The final numerical answer is extracted from the calculation.
    # total_assets = 2000000 + 4000000 + 50000 + 30000 + 20000 + 500000 = 6600000
    # total_valid_legacies = 230000
    # residuary_estate = 6600000 - 230000 = 6370000
    # "<<<6370000>>>"