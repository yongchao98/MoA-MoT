import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string stream
sys.stdout = io.StringIO()

# Step 1: Define the baseline case.
# According to the problem description, Plot 2 represents the initial, unaltered parameters.
# So, the code for plot 2 is '0'.
code_2 = '0'

# Step 2: Analyze the plots to determine the parameter change for each.

# Plot 3 Analysis:
# The top plot for case 3 shows that the concentration of reactant A (green) drops to nearly zero in the center of the slab.
# Simultaneously, the concentration of product B (purple) is also very low. This indicates that reactant A is consumed near the boundaries
# before it has a chance to diffuse into the center. This is a classic signature of a diffusion-limited process.
# This corresponds to halving the diffusion coefficient, D.
# Code: d
code_3 = 'd'

# Step 3: Rank the remaining plots (1, 4, 5, 6) by reaction speed relative to the baseline (Plot 2).
# A faster reaction will consume more reactant A. We can judge the reaction speed by looking at the minimum concentration of A (the bottom of the green curve).
# From visual inspection of the top plots:
# cA_min(Plot 2, baseline) is approx 0.25.
# cA_min(Plot 5) is the lowest (approx 0.1), so it's the fastest reaction.
# cA_min(Plot 1) is next lowest (approx 0.15), so it's the second fastest reaction.
# cA_min(Plot 4) is higher than baseline (approx 0.28), so it's a slower reaction.
# cA_min(Plot 6) is the highest (approx 0.35), so it's the slowest reaction.
# The reaction speed ranking is: Plot 5 > Plot 1 > Plot 2 (baseline) > Plot 4 > Plot 6.

# Step 4: Assign parameters based on reaction speed ranking.

# Fast Reactions (faster than baseline): Plots 1 and 5.
# The parameters that speed up the reaction are K (doubling k) and n (halving n).
# Plot 5 shows the most extreme increase in reaction speed. Doubling the rate constant K is a direct and powerful way to achieve this.
# Halving the reaction order n also speeds up the reaction but the effect is less intuitive. We assign the most extreme plot to the most direct change.
# Code for Plot 5: K
code_5 = 'K'
# Code for Plot 1: n
code_1 = 'n'

# Slow Reactions (slower than baseline): Plots 4 and 6.
# The parameters that slow down the reaction are k (halving k) and N (doubling n).
# Plot 6 shows the most extreme decrease in reaction speed. Doubling the reaction order N can have a dramatic slowing effect (since c_A < 1, c_A^N decreases significantly).
# Halving the rate constant k is a more moderate change. We assign the most extreme plot to the most dramatic change.
# Code for Plot 6: N
code_6 = 'N'
# Code for Plot 4: k
code_4 = 'k'

# Step 5: Assemble the final six-character string.
final_answer = code_1 + code_2 + code_3 + code_4 + code_5 + code_6

# This variable is used to hold the final answer, and it will be captured and displayed.
# Restore stdout and print the captured output
captured_output = sys.stdout.getvalue()
sys.stdout = stdout_backup
print(final_answer)