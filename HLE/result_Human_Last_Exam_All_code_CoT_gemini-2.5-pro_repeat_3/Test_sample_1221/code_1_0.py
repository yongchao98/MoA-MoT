def calculate_residuary_estate():
    """
    Calculates the value of Betty's residuary estate based on Hong Kong law.
    This involves:
    1. Valuing the assets in the estate.
    2. Valuing the valid legacies to be paid out.
    3. Subtracting the legacies from the assets.
    """

    # --- Step 1 & 2: Identify and Calculate the Total Value of Assets in the Estate ---
    print("--- Calculating the Total Value of Betty's Estate ---")

    # Flat A: Owned as joint tenants. Under Hong Kong's Conveyancing and Property Ordinance (Cap. 219),
    # when two people die simultaneously, the younger is deemed to have survived the elder.
    # Betty (54) is younger than Alex (57), so she is deemed to have survived him.
    # The flat passes to her by survivorship and becomes part of her estate.
    flat_a_orchid = 2000000
    print(f"Value of Flat A at Orchid Apartment (passes to Betty's estate by survivorship): {flat_a_orchid}")

    # Flat B: Owned solely by Betty.
    flat_b_peach = 4000000
    print(f"Value of Flat B at Peach Apartment (solely owned): {flat_b_peach}")

    # Other assets
    bank_accounts = 50000
    nyse_shares = 30000
    personal_items = 20000
    jewellery = 500000
    print(f"Value of cash in bank accounts: {bank_accounts}")
    print(f"Value of NYSE shares: {nyse_shares}")
    print(f"Value of personal items: {personal_items}")
    print(f"Value of jewellery collection: {jewellery}")

    # Sum of all assets in the estate
    total_assets = flat_a_orchid + flat_b_peach + bank_accounts + nyse_shares + personal_items + jewellery
    print(f"\nTotal value of assets in the estate: {flat_a_orchid} + {flat_b_peach} + {bank_accounts} + {nyse_shares} + {personal_items} + {jewellery} = {total_assets}\n")

    # --- Step 3 & 4: Analyze Legacies and Calculate Total Payouts ---
    print("--- Analyzing Legacies and Calculating Payouts ---")

    # Clause 4: Gift of Flat C. This gift is 'adeemed' because the property was sold before Betty's death.
    # The gift therefore fails.
    gift_flat_c = 0
    print(f"Payout for gift of Flat C to Cindy Tran: {gift_flat_c} (Failed due to ademption)")

    # Clause 5a: Gift to friends. This gift fails for 'uncertainty of objects' because the schedule
    # listing the friends was never created. The money falls into the residue.
    gift_friends = 0
    print(f"Payout for gift to friends: {gift_friends} (Failed due to uncertainty)")

    # Clause 5b: Gift to Wills Lawyers & Co. This is a valid pecuniary legacy.
    gift_lawyers = 230000
    print(f"Payout for gift to Wills Lawyers & Co: {gift_lawyers} (Valid)")

    # Clause 6: Gift to RSPCA. This clause was 'deleted' (revoked) by Betty. The gift is not payable.
    gift_rspca = 0
    print(f"Payout for gift to RSPCA: {gift_rspca} (Revoked)")

    # Sum of all valid legacies to be paid out
    total_payouts = gift_flat_c + gift_friends + gift_lawyers + gift_rspca
    print(f"\nTotal value of payouts from the estate: {total_payouts}\n")

    # --- Step 5: Calculate the Residuary Estate ---
    print("--- Calculating the Final Residuary Estate ---")
    residuary_estate = total_assets - total_payouts
    print(f"Residuary Estate = Total Assets - Total Payouts")
    print(f"Residuary Estate = {total_assets} - {total_payouts} = {residuary_estate}")

    return residuary_estate

if __name__ == '__main__':
    final_answer = calculate_residuary_estate()
    # The final answer is wrapped in <<<>>> as requested.
    print(f"\n<<<The total value of Betty's residuary estate is {final_answer} HKD.>>>")
