# Step 1: Determine the chromatic number of the Asian sovereign state graph before 1991.
# Before the formation of the modern Central Asian states, the graph of Asian countries, while complex,
# lacked a clear, simple subgraph that required 4 colors (a 4-chromatic subgraph like a K4 or an odd wheel).
# Therefore, we can establish its chromatic number as 3.
chi_before = 3
print(f"The chromatic number of the Asian sovereign nation graph before the Soviet dissolution, chi_before, was {chi_before}.")

# Step 2: Determine the chromatic number after 1991.
# The dissolution of the USSR created a new set of bordering countries in Central Asia.
# This formed a W6 "wheel" subgraph: Uzbekistan acts as the central hub, and its five neighbors
# (Kazakhstan, Kyrgyzstan, Tajikistan, Afghanistan, Turkmenistan) form a 5-cycle rim around it.
# A wheel graph with an odd number of rim vertices (like this W6) is not 3-colorable.
# This new structure requires 4 colors, thus increasing the chromatic number of the entire Asian map to 4.
chi_after = 4
print(f"The dissolution created a 4-chromatic subgraph, making the new chromatic number, chi_after, {chi_after}.")

# Step 3: Calculate the incremental change, delta_soviet.
# This is the difference between the chromatic number after and before the event.
delta_soviet = chi_after - chi_before
print(f"The incremental change, delta_soviet, is chi_after - chi_before = {chi_after} - {chi_before} = {delta_soviet}.")

# Step 4: Analyze the change in planarity to determine beta.
# A political map on a globe is the definition of a planar graph.
# The subdivision of one country (the USSR) into several smaller countries is an operation that preserves planarity.
# Since the planarity of the graph did not change, beta is assigned the value 1.
beta = 1
print(f"The planarity of the graph did not change, therefore beta = {beta}.")

# Step 5: Calculate and print the final requested result.
final_answer = beta * delta_soviet
print(f"The final result is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}.")