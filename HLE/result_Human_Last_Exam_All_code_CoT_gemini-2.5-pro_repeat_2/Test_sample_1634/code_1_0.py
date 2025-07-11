import sys

# We are looking for the smallest non-negative integer n such that there exists
# an n-point topological space that is not irreducible.
# A space is NOT irreducible (i.e., it is reducible) if it can be written as
# a finite union of its proper closed subsets.
# A subset Z is proper if Z is not equal to the entire space X.

print("Step 1: Analyze the case n = 0.")
print("Let X be a 0-point space, so X = {}.")
print("The only closed set is X itself.")
print("There are no proper closed subsets because there are no proper subsets at all.")
print("Therefore, a 0-point space is irreducible.")
print("-" * 40)

print("Step 2: Analyze the case n = 1.")
print("Let X be a 1-point space, e.g., X = {p}.")
print("The only proper subset of X is the empty set {}.")
print("In any topology, the empty set is always a closed set.")
print("To be reducible, X must be a union of its proper closed subsets. The only candidate is the empty set.")
print("However, the union of any number of empty sets is still the empty set, not X.")
print("Therefore, any 1-point space is irreducible.")
print("-" * 40)

print("Step 3: Analyze the case n = 2.")
print("Let X be a 2-point space. We represent the points with numbers: X = {1, 2}.")
print("To show this space can be reducible, we need to define a topology where X is a union of its proper closed subsets.")
print("Consider the discrete topology, where every subset is open. Consequently, every subset is also closed.")
print("The proper closed subsets of X are {}, {1}, and {2}.")
print("Let's choose two proper closed subsets: Z1 = {1} and Z2 = {2}.")
print("Now, we check their union to see if it equals X:")

Z1 = {1}
Z2 = {2}
X = Z1.union(Z2)

# Format the sets for printing the equation.
z1_str = str(sorted(list(Z1))).replace('[', '{').replace(']', '}')
z2_str = str(sorted(list(Z2))).replace('[', '{').replace(']', '}')
x_str = str(sorted(list(X))).replace('[', '{').replace(']', '}')

print(f"The equation demonstrating reducibility is: {z1_str} U {z2_str} = {x_str}")

print("\nSince X is the union of two proper closed subsets, this 2-point space is not irreducible.")
print("-" * 40)

print("Conclusion:")
print("For n=0 and n=1, any space is irreducible.")
print("For n=2, an example of a space that is not irreducible exists.")
final_answer = 2
print(f"The smallest non-negative integer n is {final_answer}.")

<<<2>>>