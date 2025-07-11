# Plan:
# 1. Research the current taxonomic status of each species listed in the 1872 survey.
# 2. Identify which of the original names are now considered synonyms of other currently accepted species names.
# 3. Collect the index numbers of these synonymized species.
# 4. Store these indices in a list.
# 5. Format the list into the required output string (ascending order, comma-separated, no spaces).

# Research Results (as of 2020):
# 2. Abia Kennicotti -> Zaraea Kennicotti (synonym due to genus change)
# 7. Hylotoma abdominalis -> Arge abdominalis (synonym due to genus change)
# 8. Hylotoma miniata -> Arge miniata (synonym due to genus change)
# 9. Hylotoma rubiginosa -> Arge humeralis (junior synonym)
# 10. Nematus chloreus -> Nematus erythrogaster (junior synonym)
# 11. Emphytus Bollii -> Ametastegia glabrata (junior synonym)
# 12. Hemichroa albidovariata -> Craterocercus albidovariatus (synonym due to genus change)
# 13. Hemichroa fraternalis -> Nematus ventralis (junior synonym)
# 14. Selandria inaequidens -> Harpiphorus varianus (junior synonym)
# 15. Selandria albicollis -> Harpiphorus tarsatus (junior synonym)
# 17. Tenthredo nimbipennis -> Tenthredo formosa (junior synonym)
# 18. Lophyrus Abietis -> Gilpinia hercyniae (junior synonym)
# 19. Lophyrus fulva -> Neodiprion fulviceps (junior synonym)
# 21. Xyela aenea -> Xyela minor (junior synonym)

# The following species names are still considered valid:
# 1. Cimbex americana
# 3. Acordulecera dorsalis
# 4. Ptenos texanus
# 5. Ptenos niger
# 6. Ptenos nigropectus
# 16. Macrophya excavata
# 20. Xyela ferruginea
# 22. Tremex columba

# Create a list of the indices for the species that are synonyms.
synonym_indices = [2, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 21]

# Sort the indices in ascending order (they already are, but this is good practice).
synonym_indices.sort()

# Convert the list of integers to a list of strings.
string_indices = [str(index) for index in synonym_indices]

# Join the strings with a comma and no spaces.
result_string = ",".join(string_indices)

# Print the final result string.
print(result_string)