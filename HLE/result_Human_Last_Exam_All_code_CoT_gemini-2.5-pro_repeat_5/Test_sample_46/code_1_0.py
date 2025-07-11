def explain_rmbs_crisis_root_cause():
    """
    This function explains the primary root cause behind the valuation changes
    of non-agency RMBS between 2004 and 2008.
    """

    selected_choice = "F"
    explanation_header = f"Analysis of the Correct Answer (Choice {selected_choice}): The quality of the loan issuer and RMBS originator"

    explanation_body = """
The core issue that drove the non-agency RMBS crisis was a fundamental breakdown in the quality and integrity of the loan origination and securitization process. This is the true root cause, from which other problems stemmed.

Here is the causal chain:

1.  **The Root Cause (Choice F): Poor Originator Quality:**
    The financial system shifted to an "originate-to-distribute" model. Loan issuers (like Countrywide, New Century) were incentivized to generate high volumes of mortgages because they could sell them off immediately to be packaged into RMBS. This removed their incentive to ensure the borrower could actually repay the loan. This led to a systemic decline in underwriting standards, including the creation of risky loan types (e.g., NINJA - No Income, No Job, No Assets) and lending to uncreditworthy borrowers.

2.  **The Symptoms and Enablers (Other Choices):**
    - **C. Average FICO scores:** The decline in average FICO scores was a direct symptom of the poor origination standards, not the cause itself.
    - **B. The percent of floating rate debt:** Risky adjustable-rate mortgages were a product designed by these originators to make loans appear affordable upfront, feeding the volume machine. This was a feature of the problem, not the origin.
    - **H. The rating of the credit agency:** The agencies incorrectly assigned high ratings (e.g., AAA) to these risky securities. This was a critical failure that *enabled* the sale of these RMBS to a wider market, but it did not create the underlying risk. The risk was already embedded by the originator.

3.  **The Final Trigger (The Consequence):**
    - **E. Default rates:** Ultimately, the poorly underwritten loans began to default in massive numbers, especially as teaser rates on floating-rate debt expired. High defaults were the direct mechanism that destroyed the value of the RMBS, but they were the inevitable consequence of the poor quality loans created at the origin (F).

In conclusion, while factors like default rates directly impacted the cash flows and thus the value, they were the predictable outcome of the root cause: a compromised loan origination and securitization pipeline.
"""

    print(explanation_header)
    print("=" * len(explanation_header))
    print(explanation_body)

explain_rmbs_crisis_root_cause()