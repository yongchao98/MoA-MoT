# This script determines the number of atoms transferred between compounds in a multi-step synthesis.

# Question 1: How many carbons from compound 11 are present in compound 1?
# Compound 11 (cis-2-butene-1,4-diol) has 4 carbons.
carbons_in_11 = 4
# Ozonolysis cleaves the molecule into two identical fragments. The synthesis uses one fragment.
carbons_from_11_in_use = carbons_in_11 / 2
# This 2-carbon fragment is maintained throughout the entire synthesis to compound 1.
carbons_from_11_in_1 = int(carbons_from_11_in_use)


# Question 2: How many oxygens from compound 11 are present in compound 14?
# Compound 11 has 2 oxygen atoms (from two -OH groups).
oxygens_in_11 = 2
# After the molecule is split, one fragment is used, which contains one oxygen from compound 11.
oxygens_from_11_in_use = oxygens_in_11 / 2
# This oxygen atom (in a protected -OTES group) is stable and remains present in compound 14.
oxygens_from_11_in_14 = int(oxygens_from_11_in_use)


# Question 3: How many nitrogens from compound 7 are present in compound 10?
# The synthesis order is ... -> 10 -> 7.
# Compound 10 is formed using MeNO2, giving it one nitrogen atom.
nitrogens_in_10 = 1
# This nitrogen atom is conserved when compound 10 is converted to compound 7.
# Therefore, the nitrogen atom in 10 is the same one that is later found in 7.
# The question asks for the number of N atoms in 10 that are also in 7.
nitrogens_from_7_in_10 = nitrogens_in_10


# Final result: Print the three numbers separated by commas.
print(f"The number of carbons from compound 11 in compound 1 is: {carbons_from_11_in_1}")
print(f"The number of oxygens from compound 11 in compound 14 is: {oxygens_from_11_in_14}")
print(f"The number of nitrogens from compound 7 in compound 10 is: {nitrogens_from_7_in_10}")
print(f"{carbons_from_11_in_1},{oxygens_from_11_in_14},{nitrogens_from_7_in_10}")
<<<2,1,1>>>