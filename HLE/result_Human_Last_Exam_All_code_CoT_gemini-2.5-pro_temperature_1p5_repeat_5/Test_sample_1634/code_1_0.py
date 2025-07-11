# This script determines the smallest non-negative integer n for which an n-point 
# topological space exists that is not irreducible.

# Definition: A topological space X is irreducible if it cannot be written as a finite
# union of proper closed subsets.
# We are looking for a space that is NOT irreducible (i.e., reducible), which means
# it CAN be written as such a union.

print("This script will find the smallest n for which an n-point space can be non-irreducible.\n")

print("--- Step 1: Analyze n = 0 ---")
print("Let X be the 0-point space, which is the empty set, X = {}.")
print("A 'proper' subset of X must be a subset that is not equal to X.")
print("The only subset of X is the empty set itself. Since it is equal to X, there are no proper subsets.")
print("Conclusion: The 0-point space cannot be written as a union of proper closed subsets, so it is irreducible.\n")

print("--- Step 2: Analyze n = 1 ---")
print("Let X be a 1-point space, for example, X = {p}.")
print("The only proper subset of X is the empty set, {}.")
print("The empty set is always a closed set in any topology.")
print("To be non-irreducible, X would need to be a union of its proper closed subsets.")
print("However, the union of any number of empty sets is still the empty set, which is not equal to X.")
print("Conclusion: Any 1-point space is always irreducible.\n")

print("--- Step 3: Analyze n = 2 ---")
print("Let's consider a 2-point space. We can represent the points as numbers: X = {1, 2}.")
print("To show this space can be non-irreducible, we must find a topology on X where it is a union of proper closed subsets.")
print("Let's try to express X as X = Z1 U Z2, where Z1 and Z2 are proper closed subsets.")

# We choose two proper subsets whose union is X.
Z1_val = 1
Z2_val = 2
Z1 = {Z1_val}
Z2 = {Z2_val}

print(f"Let's try to use the subsets Z1 = {Z1} and Z2 = {Z2}.")
print("For this to work, Z1 and Z2 must both be CLOSED sets in our chosen topology.")
print("A set is closed if its complement is open.")
print(f"The complement of Z1 is {{ {Z2_val} }}. This set must be open.")
print(f"The complement of Z2 is {{ {Z1_val} }}. This set must also be open.")
print("A topology where every individual point's set is open is the 'discrete topology'.")
print("In the discrete topology on X = {1, 2}, the open sets are {}, {1}, {2}, and {1, 2}.")
print("Since every subset is open, every subset is also closed.")
print(f"Therefore, Z1 = {Z1} and Z2 = {Z2} are proper closed subsets.")

print("\nNow, let's form the union of these two proper closed subsets:")
# The final equation requires printing each number
print(f"{Z1.union(Z2)} = {Z1} U {Z2}")
print("This shows that the 2-point space X can be written as a union of two proper closed subsets.")
print("Conclusion: A 2-point space can be non-irreducible.\n")


print("--- Final Conclusion ---")
print("We have shown that for n=0 and n=1, any space is irreducible.")
print("We found an example of a 2-point space that is not irreducible.")
final_answer = 2
print(f"Therefore, the smallest non-negative integer n is {final_answer}.")

<<<2>>>