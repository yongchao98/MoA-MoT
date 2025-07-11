# 1. Represent the connectivity data for the PGp area from the graph.
# The keys are the connected brain areas, and the values are their
# approximate connectivity strength, visually estimated from the polar plot (scale 0-8).
pgp_connectivity_data = {
    # Insula (light green spokes)
    'Id1': 5.5,
    'Ig2': 5.5,
    'Ig1': 5.5,
    # Parietal (red spokes)
    '7PC': 1.5,
    '7A': 1.5,
    'hIP2': 1.5,
    # Frontal (light blue spokes)
    '44': 1.5,
    'FOperc': 1.5,
    # Other areas have strengths <= 1.0 and are not listed for simplicity.
}

# 2. Define the significance threshold.
# The black circle on the plot, representing p < 0.001, is at a visual strength of 2.0.
significance_threshold = 2.0

# 3. Identify the areas with connectivity strength greater than the threshold.
# These are the "strongly connected" areas.
strongly_connected_areas = []
for area, strength in pgp_connectivity_data.items():
    if strength > significance_threshold:
        strongly_connected_areas.append(area)

# 4. Print the result. The question asks for the areas to which PGp is most
# strongly and narrowly connected. The code identifies these areas.
print("The analysis of the PGp polar plot identifies the following areas as having a connectivity strength greater than the significance threshold of {}:".format(significance_threshold))
for area in strongly_connected_areas:
    print("- " + area)

print("\nBased on this analysis, the correct option lists these three areas.")
print("The answer is the combination of: {}, {}, and {}".format(strongly_connected_areas[0], strongly_connected_areas[1], strongly_connected_areas[2]))
