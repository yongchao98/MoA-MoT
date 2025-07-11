import collections

# Step 1: Define the mapping between plot number and the integer 'n'
# This is based on ranking the plots from most stable (n=0) to most chaotic (n=8).
plot_to_n = {
    1: 0,
    3: 1,
    9: 2,
    2: 3,
    7: 4,
    8: 5,
    5: 6,
    4: 7,
    6: 8,
}

# Step 2: Calculate eccentricity 'e' for each plot number.
plot_eccentricities = {}
print("Calculating the eccentricity for each plot:")
# Sort the mapping by plot number for a clear, ordered output.
sorted_plot_to_n = collections.OrderedDict(sorted(plot_to_n.items()))

for plot_num, n in sorted_plot_to_n.items():
    eccentricity = n / 20.0
    plot_eccentricities[plot_num] = eccentricity
    print(f"Plot {plot_num}: n = {n}, eccentricity e = {n} / 20 = {eccentricity:.2f}")

# Step 3: Sort the plots by their eccentricity in ascending order.
# The `key=lambda item: item[1]` sorts the dictionary items by their value (the eccentricity).
sorted_plots = sorted(plot_eccentricities.items(), key=lambda item: item[1])

# Step 4: Extract the plot numbers from the sorted list.
ordered_plot_numbers = [item[0] for item in sorted_plots]

# Step 5: Print the final list in the required format.
final_answer_string = "{" + ", ".join(map(str, ordered_plot_numbers)) + "}"
print("\nPlot numbers in ascending order of eccentricity:")
print(final_answer_string)