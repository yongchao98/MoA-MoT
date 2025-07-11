def analyze_securities_distributions():
    """
    Analyzes several distribution scenarios against Ontario securities regulations
    as of January 2020 to determine compliance.
    """

    # --- Regulatory Rules (Simplified for this analysis) ---
    # Accredited Investor (Individual) Thresholds from NI 45-106
    ACCREDITED_INVESTOR_NET_INCOME = 200000
    ACCREDITED_INVESTOR_NET_ASSETS = 5000000
    ACCREDITED_INVESTOR_FINANCIAL_ASSETS = 1000000

    # --- Scenarios based on the user's question ---
    scenarios = [
        {
            "id": "A",
            "description": "A distribution of bonds without a prospectus by Bank of China (Canada) aimed at raising funds from retail investors.",
            "issuer_type": "bank",
            "bank_schedule": "II",
            "prospectus_used": False,
            "investor_type": "retail",
        },
        {
            "id": "B",
            "description": "A distribution of bonds by JPMorgan Chase without a prospectus aimed at raising funds from a large number of investors in Canada.",
            "issuer_type": "foreign_bank",
            "bank_schedule": None, # Not a Canadian Schedule I or II bank
            "prospectus_used": False,
            "investor_type": "large_number",
        },
        {
            "id": "C",
            "description": "A distribution of shares by a private issuer to an individual investor who has no connection to the company.",
            "issuer_type": "private_issuer",
            "prospectus_used": False,
            "investor_type": "individual",
            "investor_details": {
                "salary": 35000,
                "net_assets": 10000,
                "connection_to_company": False,
            }
        },
        {
            "id": "D",
            "description": "A distribution of shares without a prospectus by Fairstone Bank of Canada aimed at raising funds from a large number of retail investors, one of whom has zero dollars of net assets and is unemployed.",
            "issuer_type": "bank",
            "bank_schedule": "I",
            "prospectus_used": False,
            "investor_type": "retail",
        },
        {
            "id": "E",
            "description": "A distribution of shares without a prospectus by Caisse populaire acadienne ltée aimed at raising funds from a large number of retail investors, one of whom has zero dollars of net assets and is unemployed.",
            "issuer_type": "credit_union",
            "prospectus_used": False,
            "investor_type": "retail", # "Retail" implies a public offering, not one restricted to members
        },
    ]

    compliant_scenario = None

    print("Analyzing scenarios for compliance with Ontario securities regulations...\n")

    for s in scenarios:
        print(f"--- Analyzing Scenario {s['id']} ---")
        print(f"Description: {s['description']}")
        is_compliant = False
        reason = ""

        if s['prospectus_used']:
            is_compliant = True
            reason = "Compliant because a prospectus was used."
        else:
            # Check for Prospectus Exemptions
            if s['issuer_type'] == 'bank' and s['bank_schedule'] in ['I', 'II']:
                is_compliant = True
                reason = "Compliant. The 'Bank Exemption' (NI 45-106, s. 2.37) allows Canadian Schedule I and II banks to issue their own securities without a prospectus to any type of investor, regardless of their financial status."
            
            elif s['issuer_type'] == 'foreign_bank':
                is_compliant = False
                reason = "Not compliant. The issuer is not a Canadian Schedule I or II bank, so the standard Bank Exemption does not apply for a broad distribution to the public without a prospectus."

            elif s['issuer_type'] == 'private_issuer':
                investor = s['investor_details']
                # Check Accredited Investor status
                is_accredited = (investor['salary'] > ACCREDITED_INVESTOR_NET_INCOME or 
                                 investor['net_assets'] >= ACCREDITED_INVESTOR_NET_ASSETS)
                
                print("Checking Accredited Investor status for the final equation:")
                print(f"Investor Salary ${investor['salary']} vs. Threshold ${ACCREDITED_INVESTOR_NET_INCOME} -> {'Pass' if investor['salary'] > ACCREDITED_INVESTOR_NET_INCOME else 'Fail'}")
                print(f"Investor Net Assets ${investor['net_assets']} vs. Threshold ${ACCREDITED_INVESTOR_NET_ASSETS} -> {'Pass' if investor['net_assets'] >= ACCREDITED_INVESTOR_NET_ASSETS else 'Fail'}")
                
                if is_accredited:
                    is_compliant = True
                    reason = "Compliant under the Accredited Investor exemption."
                elif not investor['connection_to_company']:
                    is_compliant = False
                    reason = "Not compliant. The investor does not qualify as an 'Accredited Investor' and has no connection to the company, making the 'Private Issuer Exemption' unavailable."
            
            elif s['issuer_type'] == 'credit_union':
                is_compliant = False
                reason = "Not compliant. The 'Credit Union Exemption' (NI 45-106, s. 2.38) requires that securities are distributed *only* to members of the credit union. The term 'large number of retail investors' implies a distribution to the general public, not one restricted to members."

        print(f"Result: {'Compliant' if is_compliant else 'Not Compliant'}")
        print(f"Reasoning: {reason}\n")
        
        if is_compliant:
            compliant_scenario = s['id']

    print("--- Conclusion ---")
    print("Based on the analysis, both scenarios A and D describe compliant distributions under the Bank Exemption.")
    print("However, scenario D is a stronger example because it includes an investor with no assets, correctly demonstrating that investor financial suitability is not a requirement under this specific exemption.")
    print("Therefore, D is the best answer representing a compliant distribution.")


if __name__ == "__main__":
    analyze_securities_distributions()