def analyze_lc_documents():
    """
    Analyzes the presented documents against the Letter of Credit (LC)
    requirements to identify discrepancies based on UCP 600 rules.
    """

    print("Analyzing the Letter of Credit documents based on UCP 600 rules...\n")

    # Discrepancy Analysis
    discrepancy_a = {
        "statement": "A. The document has following discrepancy: unsigned invoice",
        "is_discrepancy": False,
        "reason": "According to UCP 600 Article 18, an invoice does not need to be signed. Therefore, an unsigned invoice is NOT a discrepancy."
    }

    discrepancy_b = {
        "statement": "B. The document has following discrepancy: Bill of lading is not presented in full set",
        "is_discrepancy": True,
        "reason": ("The Bill of Lading states that 3 originals were issued, but only 1 was presented. UCP 600 Article 20 "
                   "requires presentation of a full set of originals. This IS a discrepancy.")
    }

    discrepancy_c = {
        "statement": "C. The document has following discrepancy: Bill of lading is not made out as per LC's term",
        "is_discrepancy": True,
        "reason": ("The LC requires the B/L 'to order of issuing bank'. The B/L is consigned to 'DEF Company' (a straight B/L). "
                   "Endorsement on a straight B/L is ineffective to change the consignee as required. This IS a discrepancy.")
    }

    discrepancy_d = {
        "statement": "D. The document has following discrepancy: Packing list is presented in original instead of photocopy as per LC's term",
        "is_discrepancy": False,
        "reason": "According to UCP 600 Article 14, when a photocopy is requested, presenting an original is acceptable. This is NOT a discrepancy."
    }

    discrepancy_e = {
        "statement": "E. The document has following discrepancy: Packing list is not signed by Beneficary",
        "is_discrepancy": False,
        "reason": "The LC did not specify who must sign the packing list. Therefore, this is NOT a discrepancy."
    }
    
    # Printing the step-by-step analysis
    analysis_results = [discrepancy_a, discrepancy_b, discrepancy_c, discrepancy_d, discrepancy_e]
    
    valid_discrepancies = []
    
    for result in analysis_results:
        print(f"Checking statement: {result['statement']}")
        print(f"Result: {result['is_discrepancy']}. Reason: {result['reason']}\n")
        if result['is_discrepancy']:
            valid_discrepancies.append(result['statement'][0])
            
    print("-" * 20)
    print("Conclusion:")
    if "B" in valid_discrepancies and "C" in valid_discrepancies:
        print("Both B and C are correct discrepancies.")
        final_answer = "G"
    else:
        # This part is for logical completeness, but based on the analysis, G is correct.
        final_answer = "F" 

    print("\nThe correct statement is G, because both B and C describe valid discrepancies found in the document set.")

# Execute the analysis
if __name__ == "__main__":
    analyze_lc_documents()
    print("<<<G>>>")
