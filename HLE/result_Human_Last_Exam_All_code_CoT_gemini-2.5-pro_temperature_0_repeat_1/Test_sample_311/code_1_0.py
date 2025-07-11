# This script calculates the number of atoms transferred between compounds in a chemical synthesis.

# Question 1: How many carbons from compound 11 are present in compound 1?
# 1. Compound 11 is cis-2-butene-1,4-diol, which has 4 carbons (HO-CH2-CH=CH-CH2-OH).
# 2. Ozonolysis of the double bond cleaves the 4-carbon chain into two identical 2-carbon fragments, forming aldehyde 12.
# 3. This 2-carbon fragment is then carried through the entire synthesis (12 -> 10 -> 7 -> 13 -> 14 -> 15 -> 1).
# 4. Therefore, the final product, compound 1, contains 2 carbons that originated from compound 11.
carbons_from_11_in_1 = 2

# Question 2: How many oxygens from compound 11 are present in compound 14?
# 1. Compound 11 has 2 oxygen atoms, one in each hydroxy (-OH) group.
# 2. After protection and ozonolysis, the molecule splits. Each resulting molecule of aldehyde 12 contains one of the original oxygen atoms, now in a TESO- group.
# 3. This TESO- group is stable and is not altered in the steps leading to compound 14 (12 -> 10 -> 7 -> 13 -> 14).
# 4. Thus, compound 14 contains 1 oxygen atom that originated from compound 11.
oxygens_from_11_in_14 = 1

# Question 3: How many nitrogens from compound 7 are present in compound 10?
# 1. Compound 10 is a nitroolefin, formed using nitromethane (MeNO2). It contains 1 nitrogen atom.
# 2. Compound 7 is formed via an addition reaction where compound 10 is a reactant (10 + 6 -> 7).
# 3. This means the single nitrogen atom in compound 10 is incorporated into and becomes the single nitrogen atom in compound 7.
# 4. The question asks how many nitrogens from 7 are in 10. Since it is the same atom present in both molecules, the answer is 1.
nitrogens_from_7_in_10 = 1

# The final "equation" is the set of answers to the three questions.
# We print each number as part of the final comma-separated answer.
print(f"{carbons_from_11_in_1},{oxygens_from_11_in_14},{nitrogens_from_7_in_10}")

# The final answer is encapsulated below as requested.
# <<<2,1,1>>>