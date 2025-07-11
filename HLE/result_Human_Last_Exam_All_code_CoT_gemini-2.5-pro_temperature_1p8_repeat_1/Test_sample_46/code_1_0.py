# Define the parameters for our hypothetical non-agency RMBS
principal_balance = 1000000

# 1. Pre-Crisis Expectation
# Lenders and investors expected very low default rates, e.g., 2%
expected_default_rate = 0.02

# Calculate the expected value based on this assumption
expected_value = principal_balance * (1 - expected_default_rate)

print("Scenario 1: Pre-Crisis Valuation Assumption")
print("The calculation for the expected value is based on a low anticipated default rate.")
print(f"Equation: ${principal_balance:,.0f} * (1 - {expected_default_rate}) = ${expected_value:,.0f}")
print("-" * 50)

# 2. Crisis Reality
# During the financial crisis, default rates on subprime RMBS skyrocketed to unprecedented levels, e.g., 25% or more
crisis_default_rate = 0.25

# Calculate the new value based on the realized default rate
crisis_value = principal_balance * (1 - crisis_default_rate)

print("Scenario 2: Financial Crisis Reality")
print("The calculation for the actual value reflects the spike in defaults.")
print(f"Equation: ${principal_balance:,.0f} * (1 - {crisis_default_rate}) = ${crisis_value:,.0f}")
print("-" * 50)

# Conclusion
value_lost = expected_value - crisis_value
print(f"The change in the default rate directly caused a valuation loss of ${value_lost:,.0f}.")
print("\nThis demonstrates that the default rate was the most direct and impactful determinant of the RMBS value during the crisis.")
