def analyze_tos_clauses():
    """
    Analyzes Terms of Service clauses to identify one that is likely a contract of adhesion
    that hides a material term, violating the doctrine of reasonable expectations.
    """
    clauses = [
        {
            "id": "A",
            "description": "No-compete / no-benchmark clause.",
            "is_standard": True, # Common in B2B software
            "is_material": True,
            "is_surprising_or_hidden": False,
        },
        {
            "id": "B",
            "description": "Granting the service a license to use your content.",
            "is_standard": True, # Standard for social media/content platforms
            "is_material": True,
            "is_surprising_or_hidden": False,
        },
        {
            "id": "C",
            "description": "Late fees for non-payment.",
            "is_standard": True, # Standard for any paid subscription service
            "is_material": True,
            "is_surprising_or_hidden": False,
        },
        {
            "id": "D",
            "description": "Prohibition on scraping or automated monitoring.",
            "is_standard": True, # Extremely common clause
            "is_material": False, # Generally a reasonable restriction on use
            "is_surprising_or_hidden": False,
        },
        {
            "id": "E",
            "description": "List of prohibited acts, including a ban on researching Illinois residents.",
            "is_standard": False, # The specific Illinois clause is highly unusual
            "is_material": True, # It's a significant functional limitation
            "is_surprising_or_hidden": True, # It's buried in a long list of standard rules
        },
        {
            "id": "F",
            "description": "General abuse clause, including not using AI content to train other AIs.",
            "is_standard": True, # Becoming standard for AI services
            "is_material": True,
            "is_surprising_or_hidden": False,
        },
        {
            "id": "G",
            "description": "Permission to use your profile/actions in connection with ads.",
            "is_standard": True, # Core business model of many free social networks
            "is_material": True,
            "is_surprising_or_hidden": False,
        },
    ]

    print("Analyzing TOS clauses based on the 'reasonable expectations' doctrine...")
    print("A term is problematic if it is both 'material' (important) and 'surprising or hidden'.")
    print("-" * 60)

    best_candidate = None
    for clause in clauses:
        print(f"Clause {clause['id']}: {clause['description']}")
        print(f"  - Is it a standard or reasonably expected term? {'Yes' if clause['is_standard'] else 'No'}")
        print(f"  - Is the term surprising or hidden? {'Yes' if clause['is_surprising_or_hidden'] else 'No'}")
        
        if clause["is_material"] and clause["is_surprising_or_hidden"]:
            best_candidate = clause

    print("-" * 60)
    if best_candidate:
        print(f"Conclusion: Clause {best_candidate['id']} is the most likely candidate.")
        print("Reasoning: The ban on researching individuals in Illinois is a material limitation on the service's function.")
        print("It is hidden within a long list of otherwise standard prohibitions, making it surprising and something a user would not reasonably expect to find.")
    else:
        print("No single clause perfectly fits the criteria of being both material and hidden/surprising.")

    print("\n<<<E>>>")

analyze_tos_clauses()