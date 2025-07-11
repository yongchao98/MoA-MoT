# The following indices correspond to taxa that perform photochemical synthesis.
# 2: Gloeochaete wittrockiana (Algae - Photosynthesis)
# 3: Homo sapiens (Animal - Vitamin D synthesis)
# 5: Halapricum salinum (Archaea - Bacteriorhodopsin-based phototrophy)
# 6: Aphanothece castagnei (Cyanobacteria - Photosynthesis)
# 7: Baileya pleniradiata (Plant - Photosynthesis)
# 9: Ectothiorhodosinus mongolicus (Bacteria - Anoxygenic photosynthesis)
# 10: Chlorobaculum tepidum (Bacteria - Anoxygenic photosynthesis)
# 12: Gemmatimonas phototrophica (Bacteria - Anoxygenic photosynthesis)

# Create a list of the identified indices
qualifying_indices = [2, 3, 5, 6, 7, 9, 10, 12]

# Convert numbers to strings and join them with a comma
output = ",".join(map(str, qualifying_indices))

# Print the final result
print(output)