# The task is to identify which of the listed taxa perform photochemical synthesis.
# Based on biological classification and metabolic pathways:
# 1) Acanthella cavernosa (Sponge, Animal) - No
# 2) Gloeochaete wittrockiana (Green alga) - Yes (Oxygenic photosynthesis)
# 3) Homo sapiens (Human, Animal) - No
# 4) Riftia pachyptila (Tube worm, Animal) - No (Ignores symbiosis)
# 5) Halapricum salinum (Archaea) - Yes (Bacteriorhodopsin-based phototrophy)
# 6) Aphanothece castagnei (Cyanobacteria) - Yes (Oxygenic photosynthesis)
# 7) Baileya pleniradiata (Plant) - Yes (Oxygenic photosynthesis)
# 8) Acanthella pulchra (Sponge, Animal) - No
# 9) Ectothiorhodosinus mongolicus (Purple sulfur bacteria) - Yes (Anoxygenic photosynthesis)
# 10) Chlorobaculum tepidum (Green sulfur bacteria) - Yes (Anoxygenic photosynthesis)
# 11) Stygichthys typhlops (Cavefish, Animal) - No
# 12) Gemmatimonas phototrophica (Bacteria) - Yes (Anoxygenic photosynthesis)
# 13) Myonera garretti (Clam, Animal) - No

# The indices of the species that perform photochemical synthesis.
photochemical_synthesis_indices = [2, 5, 6, 7, 9, 10, 12]

# Format the output as a comma-separated string of the indices.
# "Remember in the final code you still need to output each number"
# The join method effectively lists each number in the final string.
result_string = ",".join(map(str, photochemical_synthesis_indices))

# Print the final result.
print(result_string)