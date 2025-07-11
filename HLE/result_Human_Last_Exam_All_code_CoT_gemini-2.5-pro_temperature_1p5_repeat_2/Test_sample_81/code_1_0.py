import textwrap

def analyze_tos_clauses():
    """
    Analyzes terms of service clauses to identify the one most likely to be a contract of adhesion
    with egregious terms or one that hides a material term of service.
    """
    clauses = {
        'A': "Customer agrees that it will not build or benchmark a competitive product or service, or copy any features, functions or graphics of the Products or Services.",
        'B': "As between you and us, you own the content and information that you submit or post to our Services, and you are only granting us and our affiliates a worldwide, transferable and sublicensable right to use, copy, modify, distribute, publish and process, information and content that you provide through our Services and the services of others, without any further consent, notice and/or compensation to you or others.",
        'C': "To the extent the jurisdiction that You are located in allows You to incur late charges for failure to pay Fees in a timely manner, any amounts arising in relation to this Agreement not paid when due will be subject to a late charge of one and one-half percent (10.5%) per month on the unpaid balance or the maximum rate allowed by law, whichever is less. Without prejudice to Your rights set out elsewhere in this Agreement, all Fees are non-refundable and payable in advance.",
        'D': "You may not use any software, device, automated process, or any similar or equivalent manual process to scrape, copy, or perform measurement, analysis, or monitoring of, any portion of the Content or Services.",
        'E': "With respect to our Products and Services, You and all Users are prohibited from engaging in the following acts: ... (v) use of any automated systems or software to extract the whole or any part of the Products and Services ... (commonly known as 'scraping'), or reverse engineer the Products and Services...",
        'F': "You must not abuse, harm, interfere with, or disrupt our services or systems — for example, by: using automated means to access content from any of our services; or using AI-generated content from our services to develop machine learning models or related AI technology.",
        'G': "You give us permission to show your username, profile picture, and information about your actions (such as likes) or relationships (such as follows) next to or in connection with accounts, ads, offers, and other sponsored content that you follow or engage with that are displayed on our Products, without any compensation to you."
    }

    print("Analyzing Terms of Service Clauses:\n")

    # While all web service TOS are effectively 'contracts of adhesion', the question asks to identify one that
    # is particularly egregious or hides a material term, violating the 'doctrine of reasonable expectations'.

    analysis = {
        'A': "Standard anti-competition clause. Common in business software agreements. Expected.",
        'B': "Standard content license grant. Necessary for a service to display user-generated content. Expected.",
        'C': "Standard late fee clause for a paid service. The rate is high but is capped by law. Expected in principle.",
        'D': "Standard anti-scraping clause to protect intellectual property. Expected.",
        'E': "A list of standard prohibitions on misuse, account sharing, and scraping. Expected.",
        'F': "A standard 'do no harm' clause, updated for the AI era. Expected.",
        'G': "A clause allowing the use of a user's identity (name, photo, actions) in advertisements without compensation. This is a material term that a reasonable user might not expect. It converts a user's social activity into a commercial endorsement. This is the most likely to violate the 'doctrine of reasonable expectations'."
    }
    
    for key, text in analysis.items():
        print(f"Clause {key}: {text}")
    
    print("\n" + "="*50)
    print("Conclusion:")
    print("The final choice is G because it describes a practice that goes beyond the reasonable expectations of a typical user.")
    print("A user clicking 'like' does not reasonably expect that this action, combined with their name and profile picture, will be used as an uncompensated endorsement in an advertisement.")
    print("Final Equation: (Analysis of Reasonable Expectations) -> Identifying a Surprising and Material Term = G")
    print("="*50)


if __name__ == "__main__":
    analyze_tos_clauses()
<<<G>>>