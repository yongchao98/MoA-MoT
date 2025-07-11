def analyze_lewis_case():
    """
    Analyzes the legal situation between Lewis and Marcel based on the Sale of Goods Act (SGA).
    """

    # --- Case Facts ---
    # The price paid by Lewis for the painting.
    price_paid = 5000

    # Key terms of the oral contract.
    contract_description = "A very detailed and large landscape painting of Algonquin Park or Hudson Bay during Autumn."
    contract_purpose = "To be the centrepiece of a living room."

    # The item that was actually delivered.
    delivered_item_description = "A very small, hastily painted picture of a nearby creek."

    # --- Legal Analysis ---
    print("Step 1: Determine if the Sale of Goods Act (SGA) applies.")
    # A contract for a custom-made item to be sold is considered a contract for the sale of goods.
    sga_applies = True
    print("Conclusion: The contract is for the sale of a 'good' (the painting), so the SGA applies.\n")

    breach_of_description = 0
    breach_of_fitness_for_purpose = 0

    print("Step 2: Check for breach of implied condition - Correspondence with Description.")
    if contract_description != delivered_item_description:
        print("Result: BREACH. The delivered painting does not match the contract's description.")
        breach_of_description = 1
    else:
        print("Result: No breach.")
    print("")

    print("Step 3: Check for breach of implied condition - Fitness for Purpose.")
    # The delivered item's attributes (small, hasty) make it unfit for the stated purpose (living room centrepiece).
    is_fit_for_purpose = False
    if not is_fit_for_purpose:
        print(f"Result: BREACH. The small painting is not fit for the stated purpose: '{contract_purpose}'.")
        breach_of_fitness_for_purpose = 1
    else:
        print("Result: No breach.")
    print("")

    print("Step 4: Conclude based on the number of breaches.")
    total_breaches = breach_of_description + breach_of_fitness_for_purpose

    # The prompt requires showing the numbers in the final equation.
    print(f"The number of major breaches of implied conditions under the SGA is: {breach_of_description} + {breach_of_fitness_for_purpose} = {total_breaches}")
    
    if total_breaches > 0:
        print(f"\nConclusion: Due to these fundamental breaches, Lewis is entitled to reject the painting")
        print(f"and recover the full purchase price of ${price_paid}.")
        print("\nTherefore, the correct answer choice is D.")
    else:
        print("\nConclusion: No breaches were found.")

# Run the analysis
analyze_lewis_case()