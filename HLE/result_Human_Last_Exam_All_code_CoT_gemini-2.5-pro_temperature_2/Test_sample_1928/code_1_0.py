# This script verifies Statement III by comparing the orbital period of 16 Cygni Bb
# with the cumulative duration of the three shortest U.S. presidential administrations.

# Define the known values in days
orbital_period_16_cygni_bb = 798.5
harrison_presidency_days = 31
garfield_presidency_days = 199
taylor_presidency_days = 492

# Sum the durations of the presidencies
total_presidency_days = harrison_presidency_days + garfield_presidency_days + taylor_presidency_days

# Print the inputs and the equation
print("Analysis of Statement III:")
print(f"Orbital period of 16 Cygni Bb: {orbital_period_16_cygni_bb} days")
print(f"Duration of William Henry Harrison's administration: {harrison_presidency_days} days")
print(f"Duration of James A. Garfield's administration: {garfield_presidency_days} days")
print(f"Duration of Zachary Taylor's administration: {taylor_presidency_days} days")

print("\nFinal Equation for Cumulative Duration:")
print(f"{harrison_presidency_days} + {garfield_presidency_days} + {taylor_presidency_days} = {total_presidency_days}")

# Print the conclusion
print(f"\nComparing the cumulative duration ({total_presidency_days} days) to the orbital period ({orbital_period_16_cygni_bb} days).")
if total_presidency_days < orbital_period_16_cygni_bb:
    print("Conclusion: The sum is less than the orbital period. Statement III is TRUE.")
else:
    print("Conclusion: The sum is not less than the orbital period. Statement III is FALSE.")
