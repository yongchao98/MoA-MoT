def analyze_legal_scenario():
    """
    This function analyzes the four legal questions based on the provided text
    and prints a step-by-step conclusion.
    """
    
    # Analysis for each of the four questions
    analysis = {
        "Question 1 (Bryan's Employment Agreement)": 
            "The non-competition and non-solicitation clauses are likely ENFORCEABLE. "
            "When restrictive covenants are part of an agreement tied to the sale of a business, "
            "courts are much more likely to enforce them to protect the company's goodwill. "
            "Given Bryan's role as CEO and former owner, the terms are likely seen as reasonable.",

        "Question 2 (Ryan's Employment Agreement)":
            "The non-competition and non-solicitation clauses are likely UNENFORCEABLE. "
            "While also tied to the sale, the reasonableness of a covenant depends on the employee's role. "
            "For a 'Shift Manager,' a province-wide non-competition clause is overly broad and not justifiable.",

        "Question 3 (New Employees' Employment Agreements)":
            "The agreements are VALID, but the non-competition clauses are UNENFORCEABLE. "
            "Under Ontario's Employment Standards Act, non-compete clauses are banned for non-executive employees. "
            "However, the unenforceable clause would be severed, leaving the rest of the contract valid.",
        
        "Question 4 (The Pickup Truck)":
            "Bryan is NOT REQUIRED to transfer the truck. His promise was a 'gratuitous promise' (a gift). "
            "For a promise to be a legally enforceable contract, there must be consideration (something of value) "
            "exchanged by both parties. Stan provided no consideration for the truck."
    }

    print("--- Legal Analysis ---")
    for question, conclusion in analysis.items():
        print(f"\n{question}:\n{conclusion}")

    final_answer = "A"
    print("\n--- Conclusion ---")
    print(f"Based on the analysis, the correct answer choice is A.")
    print(f"Final Answer: {final_answer}")

analyze_legal_scenario()