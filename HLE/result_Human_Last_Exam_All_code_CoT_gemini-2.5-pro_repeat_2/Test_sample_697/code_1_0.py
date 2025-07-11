import textwrap

def analyze_securities_regulations():
    """
    Analyzes several scenarios against Ontario securities regulations as of Jan 2020
    to determine which one is compliant.
    """
    print("Analyzing which distribution complies with Ontario securities regulations (NI 45-106 Prospectus Exemptions).")
    print("-" * 80)

    # Option A Analysis
    print("Analysis of Option A: Distribution of bonds by Bank of China (Canada).")
    analysis_a = """
    Bank of China (Canada) is a Schedule II bank under the Bank Act (Canada).
    NI 45-106, section 3.2, provides a prospectus exemption for any security issued by a bank
    listed in Schedule I, II, or III. This exemption applies to retail investors regardless of
    their financial status. This distribution appears to be compliant.
    """
    print(textwrap.indent(textwrap.dedent(analysis_a), '  '))

    # Option B Analysis
    print("Analysis of Option B: Distribution of bonds by JPMorgan Chase.")
    analysis_b = """
    The issuer 'JPMorgan Chase' is ambiguous. If it's the U.S. parent company, the Canadian
    bank exemption does not apply. If it is the Canadian Schedule III branch, the exemption
    technically applies, but a broad retail offering is not typical for a Schedule III bank's
    business model in Canada. This option is highly questionable and likely non-compliant.
    """
    print(textwrap.indent(textwrap.dedent(analysis_b), '  '))

    # Option C Analysis
    print("Analysis of Option C: Distribution by a private issuer to an unconnected individual.")
    salary = 35000
    net_assets = 10000
    accredited_income = 200000
    accredited_net_assets = 5000000
    analysis_c = f"""
    The investor's salary of ${salary:,} and net assets of ${net_assets:,} do not meet the accredited
    investor thresholds (e.g., >${accredited_income:,} income or >${accredited_net_assets:,} net assets).
    The 'private issuer' exemption requires a pre-existing connection to the company, which is
    explicitly denied in the question. This distribution is clearly non-compliant.
    """
    print(textwrap.indent(textwrap.dedent(analysis_c), '  '))

    # Option D Analysis
    print("Analysis of Option D: Distribution of shares by Fairstone Bank of Canada.")
    analysis_d = """
    Fairstone Bank of Canada is a Schedule I bank. Similar to option A, the bank exemption in
    NI 45-106 s. 3.2 applies. The exemption covers shares and applies to all retail investors,
    even one with zero net assets. This distribution appears to be compliant.
    """
    print(textwrap.indent(textwrap.dedent(analysis_d), '  '))

    # Option E Analysis
    print("Analysis of Option E: Distribution of shares by Caisse populaire acadienne ltée.")
    analysis_e = """
    A caisse populaire is a credit union. NI 45-106, section 3.5, provides an exemption for
    securities issued by a credit union, provided the distribution is 'only to members'.
    For a caisse populaire, which is a member-owned cooperative, the purchase of shares is
    the method by which an investor becomes a member. Therefore, the condition is implicitly
    satisfied. This scenario represents a fundamental and archetypal use of the exemption
    that is core to the credit union business model.
    """
    print(textwrap.indent(textwrap.dedent(analysis_e), '  '))

    # Conclusion
    print("-" * 80)
    conclusion = """
    Conclusion: Options B and C are non-compliant. Options A, D, and E all appear to be compliant
    based on the rules. However, in a 'best-fit' scenario typical of regulatory questions,
    Option E is the strongest answer. The distribution of shares to retail investors (who in
    turn become members) is the fundamental capitalization model for a caisse populaire,
    making it a classic and unambiguous example of a compliant distribution under the relevant
    exemption.
    """
    print(textwrap.dedent(conclusion))

analyze_securities_regulations()
print("<<<E>>>")