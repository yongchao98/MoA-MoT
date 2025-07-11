import math

# Define the parameters from the problem description
num_quarters = 20
quarterly_risk_increase = 0.01
agent_a_calculated_failure_p = 0.2

print("This script analyzes Agent A's reasoning error.")
print("-------------------------------------------------\n")

# --- Step 1: Replicating Agent A's Flawed Calculation ---
# Agent A's error was assuming that the total probability of failure
# is the simple sum of the risk increases for each quarter.
agent_a_calculation = num_quarters * quarterly_risk_increase
print("Agent A's flawed calculation is based on linear addition:")
print(f"{num_quarters} quarters * {quarterly_risk_increase} risk per quarter = {agent_a_calculation:.2f}\n")


# --- Step 2: Performing the Correct Calculation ---
# The correct method calculates the probability of NOT failing in any quarter (i.e., total success)
# and subtracts this from 1.
# Probability of success in one quarter = 1 - 0.01 = 0.99
# Probability of success in all 20 quarters = 0.99^20
prob_success_per_quarter = 1 - quarterly_risk_increase
total_prob_success = math.pow(prob_success_per_quarter, num_quarters)
correct_prob_failure = 1 - total_prob_success

print("The correct method based on Agent A's own premise of a 1% quarterly risk:")
print("The final equation for the cumulative probability of failure is P(F) = 1 - (P(Success))^n.")
print("Each number in the final equation is printed below:\n")

# Print the final equation with all its numeric components, as requested.
# Using format to make it clear
equation_str = "P(Failure) = 1 - (1 - {q_risk})^{n} = 1 - ({prob_s})^{n} = 1 - {total_s:.4f} = {final_p:.4f}".format(
    q_risk=quarterly_risk_increase,
    n=num_quarters,
    prob_s=prob_success_per_quarter,
    total_s=total_prob_success,
    final_p=correct_prob_failure
)
print("P(Failure) = 1 - (1 - " + str(quarterly_risk_increase) + ")^" + str(num_quarters) + " = " + str(round(correct_prob_failure, 4)))
print("\n")
print("Agent A's miscalculation (0.20) resulted from the false premise that he could linearly add the quarterly risk increases. The correct calculation, based on his own premise, would have been ~0.1821.")
print("This identifies the core flaw in his reasoning.")
