# This script solves the problem by first establishing a mapping
# between the plot grid and the simulation parameters (Da1, Da2),
# and then calculating Da3 to sort the plots.

# 1. Deduced Parameter Grid Mapping
# - Rows (top to bottom) correspond to decreasing Da_1.
# - Columns (left to right) correspond to decreasing Da_2/Da_1 ratio.
da1_map = {0: 8, 1: 4, 2: 2}  # Row index (0,1,2) to Da_1
da2_ratio_map = {0: 2.0, 1: 1.0, 2: 0.5} # Column index (0,1,2) to Da_2/Da_1 ratio

plot_data = []

print("Calculating Da_3 for each plot based on the deduced grid structure...")
print("-" * 55)

# 2. Iterate through the 3x3 grid of plots to calculate Da_3
for r in range(3):
    for c in range(3):
        plot_num = r * 3 + c + 1
        
        # Assign parameters based on the deduced mapping
        da1 = da1_map[r]
        da2 = da2_ratio_map[c] * da1
        
        # The problem defines Da_3 = Da_2 + 1/2 * Da_1
        da3 = da2 + 0.5 * da1
        
        plot_data.append({'plot': plot_num, 'da3': da3})
        
        # Output the equation for each plot as requested
        print(f"Plot {plot_num}: Da_3 = {int(da2)} + 0.5 * {int(da1)} = {int(da3)}")

# 3. Sort plots by Da_3 and generate the final answer
# Sort the plots based on the calculated Da_3 values in ascending order
sorted_plots = sorted(plot_data, key=lambda x: x['da3'])

# Extract the plot numbers in the sorted order
sorted_plot_numbers = [p['plot'] for p in sorted_plots]

# Combine the sorted plot numbers into a single 9-digit integer string
final_answer = "".join(map(str, sorted_plot_numbers))

print("-" * 55)
print("The plot numbers arranged in ascending order of Da_3 are:")
print(final_answer)