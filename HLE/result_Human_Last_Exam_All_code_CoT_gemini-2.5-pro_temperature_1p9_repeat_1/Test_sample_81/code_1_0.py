import textwrap

def analyze_tos_clauses():
    """
    Analyzes several Terms of Service clauses to identify which is most likely
    a contract of adhesion with egregious terms or one that hides material terms.
    """

    clauses = {
        'A': "Customer agrees that it will not build or benchmark a competitive product or service, or copy any features, functions or graphics of the Products or Services.",
        'B': "As between you and us, you own the content and information that you submit or post to our Services, and you are only granting us and our affiliates a worldwide, transferable and sublicensable right to use, copy, modify, distribute, publish and process, information and content that you provide through our Services...",
        'C': "To the extent the jurisdiction that You are located in allows You to incur late charges for failure to pay Fees in a timely manner, any amounts arising in relation to this Agreement not paid when due will be subject to a late charge of one and one-half percent (10.5%) per month...",
        'D': "You may not use any software, device, automated process, or any similar or equivalent manual process to scrape, copy, or perform measurement, analysis, or monitoring of, any portion of the Content or Services.",
        'E': "With respect to our Products and Services, You and all Users are prohibited from engaging in the following acts: (i) using the Products or Services for a commercial purpose; ... (vi) research or identify any individuals known to be residing in the State of Illinois, U.S.A.; and (vii) bypass security protocols...",
        'F': "You must not abuse, harm, interfere with, or disrupt our services or systems — for example, by: using automated means to access content from any of our services; or using AI-generated content from our services to develop machine learning models or related AI technology.",
        'G': "You give us permission to show your username, profile picture, and information about your actions (such as likes) or relationships (such as follows) next to or in connection with accounts, ads, offers, and other sponsored content... without any compensation to you."
    }

    analysis = {
        'A': "Standard anti-competition clause common in B2B/SaaS agreements. A user would reasonably expect this.",
        'B': "Standard content licensing grant required for any social media or content hosting service to function. Widely expected.",
        'C': "Standard late fee and payment terms for a paid service. The parenthetical number '10.5%' is likely a typo for 1.5% per month (18%/yr), but the clause defers to the legal maximum, which is standard.",
        'D': "Standard anti-scraping clause found on most websites to protect their data and resources. Expected by users.",
        'E': "This clause contains a list of mostly standard prohibitions, but it hides a highly specific and material term. A user would not reasonably expect the prohibition found in item (vi), which forbids research on individuals in a single specific state (Illinois). This is likely an attempt to avoid liability under a specific state law (like BIPA) and is a material, unexpected limitation on the service.",
        'F': "A modern, but now standard, anti-scraping and anti-AI-training clause. Increasingly expected by users as of 2024.",
        'G': "The fundamental model for most free social media services (like Facebook/Instagram). This describes how they use data for ads, and it is a widely understood and expected practice."
    }
    
    print("--- Analysis of Terms of Service Clauses ---")
    wrapper = textwrap.TextWrapper(width=80, initial_indent="    ", subsequent_indent="    ")
    
    for option, text in clauses.items():
        print(f"\n[Option {option}]")
        print(wrapper.fill(text))
        print("\n  Analysis:")
        print(f"    {analysis[option]}")

    print("\n--- Conclusion ---")
    print("While all online TOS are technically contracts of adhesion, Option E stands out.")
    print("It hides a material term that violates the doctrine of reasonable expectations.")
    print("A user would not expect a geographic restriction against researching people in Illinois.")
    print("The key part of the clause is numbered (vi), which is an unusual term, hidden in a list of otherwise standard prohibitions like item (vii).")
    print("\nFinal Answer: E")

if __name__ == '__main__':
    analyze_tos_clauses()