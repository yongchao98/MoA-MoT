# The user wants to find the best corporate structure.
# Let's analyze the requirements and the options.

# Requirements:
# 1. Equal control for Alex and Tyler.
# 2. Alex to be paid via salary.
# 3. Tyler to be paid via dividends.
# 4. Option for future investors with no voting rights.

# Analysis of Options:
# Option A:
# 1. Equal Control: Yes. Sole two directors and equal voting shares (50 vs 50).
# 2. Alex Salary: Yes. As a director/employee, he can receive a salary.
# 3. Tyler Dividends: Yes. His Class A shares are dividend-eligible.
# 4. Non-Voting Investors: Yes. Class C shares are authorized and are non-voting. This provides the structural option.
# Verdict: Meets all requirements.

# Option B:
# 3. Tyler Dividends: No. Tyler gets Class A shares which are not eligible for dividends.
# Verdict: Fails a key requirement.

# Option C:
# This option is less ideal than A. It gives Alex non-dividend paying shares, which is very restrictive for a founder, even if he currently prefers a salary. While it technically meets the requirements, it's a poorly designed structure for the founders. Option A is superior as it treats both founders equitably.

# Option D:
# 4. Non-Voting Investors: No. It does not authorize any class of non-voting shares.
# Verdict: Fails a key requirement.

# Option E:
# 1. Equal Control: No. The board is 4 people: Alex (1 vote) vs. Tyler and his two sons (3 votes). This is not equal control over operations.
# Verdict: Fails a key requirement.

# Conclusion: Option A is the only one that satisfies all the requirements in a logical and equitable manner.

correct_answer = "A"
explanation = "Option A is the correct choice because it meets all four requirements stated by the business partners. \n"\
"1) Equal Control: Alex and Tyler are the sole directors and control an equal number of voting shares, ensuring equal control. \n"\
"2) Alex's Salary: The structure allows Alex to be a director and employee, enabling him to receive a salary. \n"\
"3) Tyler's Dividends: Tyler's Class A shares are eligible for dividends, satisfying his payment preference. \n"\
"4) Future Investors: The authorization of non-voting Class C shares provides the option to bring in future investors without giving them control over the company's operations."

print(f"The correct option is: {correct_answer}")
print("\nExplanation:")
print(explanation)
print("\nTherefore, the final answer is <<<A>>>")