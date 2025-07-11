import sys
from io import StringIO

# A helper function to format the output
def print_step(step_num, explanation, equation_str=None, result=None):
    """Prints a step in the derivation."""
    print(f"Step {step_num}: {explanation}")
    if equation_str:
        if result is not None:
            # We are asked to output each number in the final equation.
            # Here, we will print the components of the equation.
            parts = equation_str.split("=")
            vars_and_vals = parts[0].strip()
            calc = parts[1].strip()
            print(f"  Equation: {vars_and_vals} = {calc} = {result}")
        else:
            print(f"  {equation_str}")
    print("-" * 20)

# Main explanation
print("This problem asks for the limit of a probability in a branching random walk as a parameter h -> 0.")
print("The parameter h controls both the environment's randomness and the branching rate.")
print("The solution can be found by analyzing the behavior of the system in the limit.")
print("-" * 20)

# Step 1: Limiting Environment
print_step(1, "As h -> 0, the probability of any site being 'red' is h, which goes to 0. "
           "The environment almost surely becomes uniformly 'blue'.")

# Step 2: Limiting Random Walk
print_step(2, "In an all-blue environment, the walk has fixed jump probabilities.")
p_L_blue = "1/5"
p_R_blue = "4/5"
print_step(2.1, f"  Probability of jumping left from a blue site: p_L = {p_L_blue}")
print_step(2.2, f"  Probability of jumping right from a blue site: p_R = {p_R_blue}")

# Step 3: Drift Calculation
p_R_val = 4/5
p_L_val = 1/5
drift_v = p_R_val - p_L_val
print_step(3, "The drift of a single particle in this environment is calculated.",
           equation_str=f"v = p_R - p_L = {p_R_blue} - {p_L_blue}", result=drift_v)
print("  A positive drift (v > 0) means particles are, on average, pulled to the right (towards +infinity).")
print("-" * 20)


# Step 4: Behavior of the Branching Process
print_step(4, "For a branching random walk, if the underlying walk has a positive drift, "
           "the entire cloud of particles drifts to the right.")
print_step(4.1,"  The minimal position of the particle cloud (m_n) also tends to +infinity.")

# Step 5: Final Conclusion
print_step(5, "If the entire particle population drifts to +infinity, any fixed site (like site 0) "
           "can only be visited a finite number of times.")
print_step(6, "Therefore, the probability of site 0 being visited by infinitely many different particles is 0.")

# Final Answer formatting
# Use a string buffer to capture the print output for the final answer format
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

# "Final Equation" requested by the prompt
print("Let P(E) be the probability of the event.")
# The "final equation" is the limit itself. We state the value.
print("lim_{h->0} P(E) = 0")
sys.stdout = old_stdout
result_str = mystdout.getvalue()
# print(result_str) # This would print the captured output

# Final Answer as requested format
final_answer = 0
print("\nFinal calculation leads to the result:")
print(f"{final_answer}")