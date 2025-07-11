# Agent A's task parameters
total_duration_years = 5
quarters_per_year = 4
prob_increase_per_quarter = 0.01

# Calculate the total number of quarters
num_quarters = total_duration_years * quarters_per_year

# Agent A's flawed reasoning: He calculates the final probability by linearly adding the 
# assumed probability increase for each quarter. This is the core of his miscalculation.
# This calculation directly leads to his conclusion of p = 0.2.
flawed_final_probability = num_quarters * prob_increase_per_quarter

print("Agent A's flawed calculation was based on the false premise that cumulative probability can be found by simply adding the risk increases from each period.")
print("The calculation performed was:")
# The format f'{num_quarters:g}' ensures that the number is printed without trailing zeros if it's a whole number.
print(f"{num_quarters:g} quarters * {prob_increase_per_quarter:g} increase/quarter = {flawed_final_probability:g}")
print("\nThis method is incorrect because probabilities of failure do not add up in this linear fashion.")
