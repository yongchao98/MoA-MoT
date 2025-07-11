# Number of possible values for the first invariant, based on (x+y) mod 3.
# The invariant is a pair of binary values ((n0-n1)%2, (n1-n2)%2).
num_classes_invariant_1 = 2 * 2

# Number of possible values for the second invariant, based on (x-y) mod 3.
# The invariant is a pair of binary values ((m0-m1)%2, (m1-m2)%2).
num_classes_invariant_2 = 2 * 2

# The total number of equivalence classes is the product of the number of
# possible values for each independent invariant.
total_classes = num_classes_invariant_1 * num_classes_invariant_2

print(f"The number of possible values for the first invariant is {num_classes_invariant_1}.")
print(f"The number of possible values for the second invariant is {num_classes_invariant_2}.")
print(f"The total number of equivalence classes is the product of these two numbers.")
print(f"{num_classes_invariant_1} * {num_classes_invariant_2} = {total_classes}")
