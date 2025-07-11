# 2024 Tax Year Contribution Limits for someone under age 50
IRA_LIMIT = 7000
K401_EMPLOYEE_LIMIT = 23000
K401_TOTAL_LIMIT = 69000  # Employee + Employer + After-tax contributions

# User's current contributions
user_ira_contribution = 5000
user_k401_employee_contribution = 23000

# User's employer match details
employer_match_rate = 0.50
# Note: The prompt says "up to 50% on my 401k". We assume this means 50% of the employee's contribution, which is a common setup.
employer_match_contribution = employer_match_rate * user_k401_employee_contribution

# --- Calculations ---

# 1. Calculate remaining IRA contribution room
remaining_ira_contribution = IRA_LIMIT - user_ira_contribution

# 2. Check if standard 401(k) employee contributions are maxed out
# This will be 0, as the user has already contributed the maximum.
remaining_k401_employee_contribution = K401_EMPLOYEE_LIMIT - user_k401_employee_contribution

# 3. Calculate remaining room for after-tax 401(k) contributions (Mega Backdoor Roth)
# This is only possible if the user's 401(k) plan allows for after-tax contributions.
total_current_k401_contributions = user_k401_employee_contribution + employer_match_contribution
remaining_after_tax_k401_room = K401_TOTAL_LIMIT - total_current_k401_contributions

# 4. Calculate total potential additional contributions
total_potential_contribution = remaining_ira_contribution + remaining_after_tax_k401_room

# --- Output ---
print("Based on 2024 tax law for individuals under age 50:")
print("-" * 50)
print(f"Your standard employee 401(k) contribution is maxed out at ${int(user_k401_employee_contribution)}.")
print(f"Your employer match adds an additional ${int(employer_match_contribution)}.")
print("-" * 50)
print("Here is the breakdown of how much more you can potentially contribute:\n")

print(f"1. Remaining IRA Contribution: ${int(remaining_ira_contribution)}")
print(f"2. Remaining After-Tax 401(k) Contribution: ${int(remaining_after_tax_k401_room)}")
print("   (Note: This is only possible if your 401(k) plan specifically allows after-tax contributions. Please check with your plan administrator.)\n")

print("Final Potential Contribution Equation:")
# Final output showing each number in the equation
print(f"${int(remaining_ira_contribution)} (IRA) + ${int(remaining_after_tax_k401_room)} (After-Tax 401(k)) = ${int(total_potential_contribution)}")
