import textwrap

def analyze_case():
    """
    Analyzes the legal dispute between Lewis and Marcel
    regarding the painting and the Sale of Goods Act (SGA).
    """

    # --- Case Facts ---
    buyer = "Lewis"
    seller = "Marcel"
    price = 5000
    delivery_date = "March 1, 2022"

    contract_specifications = {
        "description": "A very detailed and large landscape painting of Algonquin Park or Hudson Bay during Autumn.",
        "purpose": "To hang as the centrepiece of their living room."
    }

    delivered_good = {
        "description": "A very small, hastily painted picture that clearly resembled a nearby creek.",
        "fitness_for_purpose": False
    }

    # --- Legal Analysis ---
    print("Legal Analysis of Lewis v. Marcel")
    print("=" * 40)

    # Step 1: Determine if the contract is for goods or services.
    print("\n1. Applicability of the Sale of Goods Act (SGA):")
    analysis1 = """\
    The contract is for a custom-made painting. While this involves the artist's skill (a service), the primary purpose of the contract is the transfer of ownership of a tangible item (the painting, a 'good'). In such hybrid contracts, courts look at the substance. Here, the substance is the sale of a good. Therefore, the SGA and its implied conditions apply. This makes option B incorrect."""
    print(textwrap.fill(analysis1, width=80))

    # Step 2: Analyze the 'fitness for purpose' condition.
    print("\n2. Breach of Implied Condition of Fitness for Purpose:")
    analysis2 = f"""\
    The SGA implies a condition that goods are fit for a particular purpose if:
    a) The buyer makes the purpose known. Lewis clearly stated he wanted a '{contract_specifications['description']}' for the specific purpose of being a '{contract_specifications['purpose']}'.
    b) The buyer relies on the seller's skill. Lewis relied on Marcel, a 'highly regarded artist'.
    c) The seller deals in goods of that description. Marcel is an artist who sells paintings.
    
    The delivered painting ({delivered_good['description']}) did not match the specified purpose. Therefore, the implied condition was breached. This makes options A and C incorrect."""
    print(textwrap.fill(analysis2, width=80))

    # Step 3: Address the nature of the oral agreement.
    print("\n3. Oral Agreement and Implied Terms:")
    analysis3 = """\
    The conditions and warranties in the SGA are 'implied' by law, meaning they are automatically part of the contract unless the parties expressly agree to exclude them. The fact that the agreement was oral and did not mention the SGA is irrelevant. The terms are included by default. This makes option E incorrect."""
    print(textwrap.fill(analysis3, width=80))

    # Step 4: Determine the remedy.
    print("\n4. Remedy for Breach of Condition:")
    analysis4 = f"""\
    A 'condition' is a fundamental term of the contract. A breach of condition gives the innocent party (Lewis) the right to repudiate (reject) the contract and claim a full refund. Lewis is entitled to the return of the ${price} he paid.
    
    The final equation for Lewis's recovery is:
    Amount Paid - Value of Non-Conforming Good Accepted = Refund
    Since Lewis can reject the good entirely, the equation is:
    ${price} - $0 = ${price}
    """
    print(textwrap.fill(analysis4, width=80))
    
    # --- Conclusion ---
    print("\n" + "=" * 40)
    print("Conclusion:")
    conclusion = """\
    A lawyer would advise Lewis that his interpretation is correct. The SGA applies, the implied condition of fitness for purpose was breached, and he is entitled to a full refund of the money he paid. Option D is the correct statement of the law."""
    print(textwrap.fill(conclusion, width=80))
    print("=" * 40)


if __name__ == "__main__":
    analyze_case()