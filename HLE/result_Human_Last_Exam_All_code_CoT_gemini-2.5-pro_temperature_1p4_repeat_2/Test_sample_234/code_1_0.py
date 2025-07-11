# Let's represent the truth value of each property with 1 for true and 0 for false.
# Based on the analysis:
# - Open: Always true
# - Closed: Not always true
# - Connected: Not always true
# - Compact: Not always true
# - Dense: Always true
# - Connected complement: Not always true
# - Trivial first singular homology group: Always true

p_open = 1
p_closed = 0
p_connected = 0
p_compact = 0
p_dense = 1
p_connected_complement = 0
p_trivial_h1 = 1

# Calculate the total number of properties that must always be true
total_true_properties = p_open + p_closed + p_connected + p_compact + p_dense + p_connected_complement + p_trivial_h1

print(f"{p_open} (Open) + {p_closed} (Closed) + {p_connected} (Connected) + {p_compact} (Compact) + {p_dense} (Dense) + {p_connected_complement} (Connected complement) + {p_trivial_h1} (Trivial H1) = {total_true_properties}")
