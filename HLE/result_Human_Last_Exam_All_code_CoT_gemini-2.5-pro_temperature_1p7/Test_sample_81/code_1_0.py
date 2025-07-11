import textwrap

def analyze_contract_clauses():
    """
    Analyzes contract clauses to identify the most likely example of a
    contract of adhesion or a clause that hides material terms.
    """
    clauses = {
        'A': "Customer agrees that it will not build or benchmark a competitive product or service...",
        'B': "As between you and us, you own the content... you are only granting us... a worldwide, transferable and sublicensable right to use, copy, modify...",
        'C': "any amounts... not paid when due will be subject to a late charge of one and one-half percent (10.5%) per month...", # Note: The prompt has a typo; 1.5% is not 10.5%. Analysis assumes a standard 1.5%.
        'D': "You may not use any software, device, automated process... to scrape, copy, or perform measurement...",
        'E': "prohibited from engaging in... (vi) research or identify any individuals known to be residing in the State of Illinois, U.S.A....",
        'F': "You must not abuse... for example, by... using AI-generated content from our services to develop machine learning models...",
        'G': "You give us permission to show your username, profile picture, and information about your actions (such as likes)... next to or in connection with... ads... without any compensation to you."
    }

    # Scoring based on 'unexpectedness' and 'egregiousness' (scale of 1-10)
    # A high score indicates a term a user would not reasonably expect or a term that is highly unfair.
    scores = {
        'A': {'unexpectedness': 2, 'egregiousness': 3},  # Standard in B2B
        'B': {'unexpectedness': 3, 'egregiousness': 4},  # Standard for content platforms
        'C': {'unexpectedness': 1, 'egregiousness': 1},  # Standard late fee clause
        'D': {'unexpectedness': 1, 'egregiousness': 2},  # Standard anti-scraping
        'E': {'unexpectedness': 9, 'egregiousness': 7},  # Highly unusual and specific legal burden shifting
        'F': {'unexpectedness': 4, 'egregiousness': 3},  # Increasingly standard anti-AI training
        'G': {'unexpectedness': 8, 'egregiousness': 9}   # Uses user identity for ads with no compensation
    }

    best_choice = None
    max_score = -1

    print("Analyzing contract clauses...\n")

    for choice, score_data in scores.items():
        unexpectedness = score_data['unexpectedness']
        egregiousness = score_data['egregiousness']
        total_score = unexpectedness + egregiousness
        if total_score > max_score:
            max_score = total_score
            best_choice = choice

    print(f"Analysis complete. The most problematic clause is '{best_choice}'.\n")
    print("Reasoning:")
    
    explanation = (
        f"Clause {best_choice} is the most likely answer. It permits the service to use a person's identity "
        "(username, profile picture) and their actions (likes, follows) as endorsements in advertisements, "
        "all without providing any compensation. This constitutes a significant appropriation of personal "
        "identity for commercial gain. A reasonable person would not expect their personal social "
        "interactions to be turned into free advertising for the platform's paying customers. "
        "This type of term exploits the unequal bargaining power inherent in a 'take-it-or-leave-it' "
        "adhesion contract and violates the doctrine of reasonable expectations."
    )
    
    print(textwrap.fill(explanation, width=80))
    
    # Fulfilling the requirement to show the final equation
    final_unexpectedness = scores[best_choice]['unexpectedness']
    final_egregiousness = scores[best_choice]['egregiousness']
    
    print("\nScoring Equation for the winning choice:")
    print(f"Unexpectedness Score: {final_unexpectedness}")
    print(f"Egregiousness Score: {final_egregiousness}")
    print(f"Final Equation: {final_unexpectedness} + {final_egregiousness} = {max_score}")
    
    print(f"\nThe answer is: {best_choice}")

analyze_contract_clauses()