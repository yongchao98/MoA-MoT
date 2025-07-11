import sys
# Make sure we are using Python 3.8+ for the walrus operator, otherwise use standard assignment
# This is just for cleaner printing format. It can be written without it.
if sys.version_info < (3, 8):
    print("This script is best viewed with Python 3.8+ for formatted printing.")

# Given parameters
d = 5
e1 = 3
e2 = 2
q = 4

# The proportion of irreducible (e1, e2)-stingray duos is given by the formula:
# P_irr = (1 - 1/q) * (1 - P_reducible_by_subspace_containment)
# where P_reducible_by_subspace_containment = P(U2=F1) + P(U1=F2) - P(U2=F1 and U1=F2)
# P(U2=F1) = 1/q^(e1*e2)
# P(U1=F2) = 1/q^(e1*e2)
# P(U2=F1 and U1=F2) = 1/q^(2*e1*e2)

e1e2 = e1 * e2

# Probability of U2 = F1
p_U2_eq_F1 = 1 / q**e1e2
# Probability of U1 = F2
p_U1_eq_F2 = 1 / q**e1e2
# Probability of both
p_both = 1 / q**(2 * e1e2)

# Probability of reducible due to U1=F2 or U2=F1
p_red_contain = p_U2_eq_F1 + p_U1_eq_F2 - p_both

# The final formula for the proportion of irreducible duos
proportion = (1 - 1/q) * (1 - p_red_contain)

# Output the explanation and calculation
print("The final answer is composed of three parts: (a), (b), and (c).")
print("(a) No, not all (3, 2)-stingray duos are irreducible.")
print("(b) Reducibility is caused by any of the conditions { (1) F1 cap F2 != {0}, (2) U1 = F2, (3) U2 = F1 } being true.")
print("(c) The proportion of irreducible (3,2)-stingray duos is calculated as follows:\n")

print(f"Let q = {q}, e1 = {e1}, e2 = {e2}. The product e1*e2 = {e1e2}.")
print("The probability of reducibility due to subspace equality is P(E) = P(U2=F1) + P(U1=F2) - P(both).")
print(f"P(U2=F1) = 1/q^(e1*e2) = 1/{q**e1e2}")
print(f"P(U1=F2) = 1/q^(e1*e2) = 1/{q**e1e2}")
print(f"P(E) = 2/q^6 - 1/q^12 = 2/{q**6} - 1/{q**12} = {p_red_contain}")
print("\nThe proportion of irreducible duos is P_irr = (1 - 1/q) * (1 - P(E)).")
print(f"P_irr = (1 - 1/{q}) * (1 - {p_red_contain})")
# Using := requires Python 3.8+
try:
    eval("term1 = 1 - 1/q")
    eval("term2 = 1 - p_red_contain")
    print(f"P_irr = ({term1}) * ({term2})")
except SyntaxError:
    term1 = 1 - 1/q
    term2 = 1 - p_red_contain
    print(f"P_irr = ({term1}) * ({term2})")

print(f"P_irr = {proportion}")
print(f"As a fraction, this is {proportion.as_integer_ratio()[0]} / {proportion.as_integer_ratio()[1]}")
final_answer_c = f"{proportion:.7f}"
print(f"\nFinal answer string: <<<(a) No (b) {{ (1), (2), (3) }} (c) {final_answer_c}>>>")
