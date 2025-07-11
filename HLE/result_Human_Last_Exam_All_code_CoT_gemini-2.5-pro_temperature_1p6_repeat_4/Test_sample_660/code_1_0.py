import sys
import io

# Redirect stdout to capture print output for the final answer format
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

print("This script determines the set of values for theta for which more truthful companies choose more lenient auditors.")
print("The logic proceeds in five steps.\n")

print("Step 1: Define the Firm's Expected Cost")
print("The expected cost C(theta, x) for a firm of type theta choosing auditor x is the sum of its tax liability and any expected penalty.")
print("The probability of an auditor refusing a report is x * (1 - theta).")
print("If the auditor refuses, the Revenue Service audits with probability p(x), and the firm pays a penalty F(theta).")
print("The cost function is: C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta).")
print("Since refusal only happens for x=1, the only relevant audit probability is p(1).")
print("Thus, C(theta, x) = t(theta) + x * (1 - theta) * p(1) * F(theta).\n")


print("Step 2: Analyze the Auditor Choice for Each Firm Type")
print("A firm's choice is determined by comparing C(theta, 0) with C(theta, 1).\n")

print("--- Analysis for a truthful firm (theta = 1):")
print("Cost with lenient auditor (x=0): C(1, 0) = t(1) + 0 * (1-1) * p(1) * F(1)")
print("C(1, 0) = t(1)")
print("Cost with strict auditor (x=1): C(1, 1) = t(1) + 1 * (1-1) * p(1) * F(1)")
print("C(1, 1) = t(1)")
print("Result: C(1, 0) is equal to C(1, 1). The truthful firm is indifferent.")
print("For this analysis, we assume a tie-breaking rule: if indifferent, the firm chooses the lenient auditor (x=0). So, the choice is x*(1) = 0.\n")

print("--- Analysis for a firm with malpractice (theta = 0):")
print("Cost with lenient auditor (x=0): C(0, 0) = t(0) + 0 * (1-0) * p(1) * F(0)")
print("C(0, 0) = t(0)")
print("Cost with strict auditor (x=1): C(0, 1) = t(0) + 1 * (1-0) * p(1) * F(0)")
print("C(0, 1) = t(0) + p(1) * F(0)")
print("Result: The choice depends on the sign of the term p(1) * F(0).")
print("The firm chooses the strict auditor (x=1) if C(0, 1) < C(0, 0), which requires:")
print("t(0) + p(1) * F(0) < t(0)  ==>  p(1) * F(0) < 0.\n")

print("Step 3: Formalize the Condition in the Question")
print("The question is: 'What is the set of values theta for which companies keeping more truthful accounts choose more lenient auditors?'")
print("This implies a comparison: truthful firms (theta=1) choose leniently (x=0) while malpractice firms (theta=0) choose strictly (x=1).")
print("So we need to find when x*(1) = 0 AND x*(0) = 1.\n")

print("Step 4: Derive Parameter Constraints")
print("From our analysis in Step 2:")
print("- x*(1) = 0 holds due to our tie-breaking assumption for the indifferent firm.")
print("- x*(0) = 1 holds if and only if p(1) * F(0) < 0.")
print("Since p(1) is a probability, it must be non-negative. For the product to be negative, we need p(1) > 0 and F(0) < 0.")
print("This means the penalty F(0) for malpractice must actually be a negative value (a subsidy), and the probability of a secondary audit p(1) must be positive.\n")

print("Step 5: Identify the Final Set for theta")
print("When the condition p(1) * F(0) < 0 is met:")
print("- The choice for a firm with theta=1 is the lenient auditor (x=0).")
print("- The choice for a firm with theta=0 is the strict auditor (x=1).")
print("The question asks for the set of theta values for which the firm exhibits the 'more lenient choice' behavior.")
print("In this specific scenario, the firm with theta=1 is the one choosing the lenient auditor, while the firm with theta=0 is not.")
print("Therefore, the set of values of theta that satisfies the condition is {1}.")

# Extract the final number for the answer equation
final_answer_set = {1}
theta_value = list(final_answer_set)[0]

sys.stdout = old_stdout
# Get the explanation text
explanation_text = captured_output.getvalue()
print(explanation_text)

print("\nFinal Answer Equation:")
print("Let S be the set of values of theta that satisfy the condition in the question.")
print(f"S = {{{theta_value}}}")

# Final Answer Block
print(f"\n<<<{theta_value}>>>")