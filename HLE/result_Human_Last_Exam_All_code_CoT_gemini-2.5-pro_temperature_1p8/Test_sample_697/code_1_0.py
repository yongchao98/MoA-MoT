def analyze_securities_distributions():
    """
    Analyzes several securities distribution scenarios for compliance
    with Ontario securities regulations as of January 2020.
    """

    # --- Accredited Investor Thresholds for an Individual ---
    # As per NI 45-106 Prospectus Exemptions
    accredited_investor_net_financial_assets = 1_000_000
    accredited_investor_net_assets = 5_000_000
    accredited_investor_income = 200_000

    scenarios = {
        "A": {
            "issuer_type": "Schedule III Bank (Bank of China (Canada))",
            "security_type": "Bonds",
            "investor_profile": "Retail investors",
            "prospectus": False
        },
        "B": {
            "issuer_type": "Foreign Parent Company (JPMorgan Chase)",
            "security_type": "Bonds",
            "investor_profile": "Large number of investors",
            "prospectus": False
        },
        "C": {
            "issuer_type": "Private Issuer",
            "security_type": "Shares",
            "investor_profile": {
                "salary": 35_000,
                "net_assets": 10_000,
                "connection": "none"
            },
            "prospectus": False
        },
        "D": {
            "issuer_type": "Schedule I Bank (Fairstone Bank of Canada)",
            "security_type": "Shares",
            "investor_profile": "Retail, including one with zero net assets and unemployed",
            "prospectus": False
        },
        "E": {
            "issuer_type": "Out-of-Province Credit Union (Caisse populaire acadienne ltée)",
            "security_type": "Shares",
            "investor_profile": "Retail, including one with zero net assets and unemployed",
            "prospectus": False
        }
    }

    print("--- Analysis of Securities Distribution Scenarios in Ontario ---")

    # --- Analysis of Scenario C ---
    c_data = scenarios["C"]["investor_profile"]
    is_accredited_C = (c_data["salary"] >= accredited_investor_income or
                       c_data["net_assets"] >= accredited_investor_net_assets)
    print("\n[Analysis of C]: Distribution by a private issuer.")
    print(f"The individual investor has a salary of ${c_data['salary']} and net assets of ${c_data['net_assets']}.")
    print(f"The income threshold for an accredited investor is over ${accredited_investor_income}, and the net asset threshold is ${accredited_investor_net_assets}.")
    print(f"The investor has no connection to the company and does not meet the financial thresholds to be an accredited investor.")
    print("Result for C: This distribution is NOT compliant, as it does not qualify for the private issuer exemption.")

    # --- Analysis of Scenario B ---
    print("\n[Analysis of B]: Distribution by JPMorgan Chase.")
    print("The prospectus exemption for banks applies to banks listed in Schedule I, II, or III of Canada's Bank Act.")
    print("JPMorgan Chase, as a U.S. parent company, is not a bank on these schedules. The exemption does not apply to it.")
    print("Result for B: This distribution is NOT compliant without a prospectus.")

    # --- Analysis of Scenarios A, D, E ---
    print("\n[Analysis of A, D, E]: These involve issuers eligible for prospectus exemptions.")
    print("The Ontario Securities Act provides a prospectus exemption for securities issued by banks (Schedule I, II, or III) and credit unions from any Canadian jurisdiction.")

    # --- Deeper look at A, D, E ---
    print("\n[Analysis of A]: Bank of China (Canada) is a Schedule III bank.")
    print("Its distribution of bonds is exempt from the prospectus requirement. This scenario appears compliant.")

    print("\n[Analysis of D]: Fairstone Bank of Canada is a Schedule I bank.")
    print("Its distribution of shares is exempt. However, selling shares (a risk asset) to a person with zero assets for fundraising raises significant suitability and public interest concerns, making the overall practice potentially non-compliant.")
    print("Result for D: This distribution is Questionable.")
    
    print("\n[Analysis of E]: Caisse populaire acadienne ltée is a credit union.")
    print("Its distribution of shares is exempt. A credit union ('caisse populaire') is a member-owned cooperative.")
    print("Their 'shares' are often nominal-value membership shares required to join and access services.")
    print("Selling a membership share to a person with zero net assets is a normal, standard, and socially valuable practice for a credit union. It is not an unsuitable investment sale.")
    print("Result for E: This distribution is FULLY compliant with the letter and spirit of the regulations.")
    
    print("\n--- Final Conclusion ---")
    print("Comparing the compliant or potentially compliant options (A, D, E), option E describes a scenario that is not only technically exempt but represents standard, compliant operating procedure for a credit union. Therefore, it is the best answer.")


if __name__ == '__main__':
    analyze_securities_distributions()
    print("<<<E>>>")