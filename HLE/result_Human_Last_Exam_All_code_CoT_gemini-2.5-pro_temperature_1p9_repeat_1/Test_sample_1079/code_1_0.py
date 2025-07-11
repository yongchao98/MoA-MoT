# The problem is to find the number of elements in the H3 reflection group
# that have a regular eigenvector with an eigenvalue of order 10.

# 1. Identify the properties of the H3 group.
# The degrees of the fundamental invariants for H3 are 2, 6, and 10.
degrees = [2, 6, 10]

# The order of the H3 group, |H3|, is the product of its degrees.
group_order = 1
for d in degrees:
    group_order *= d

# The Coxeter number, h, for H3 is its largest degree.
coxeter_number = max(degrees)

# 2. Apply the relevant mathematical theorem.
# A theorem by Springer states that an element 'g' in an irreducible Coxeter group
# has an eigenvalue with order equal to the Coxeter number 'h' if and only if
# 'g' is conjugate to a Coxeter element 'c'.
# This means we need to find the size of the conjugacy class of the Coxeter elements.

# 3. Calculate the size of the conjugacy class.
# The size of this class is given by the formula: |H3| / |C(c)|,
# where |C(c)| is the order of the centralizer of a Coxeter element.
# The order of the centralizer of a Coxeter element is the Coxeter number, h.
centralizer_order = coxeter_number

# The number of such elements is the group order divided by the centralizer order.
num_elements = group_order // centralizer_order

# 4. Print the final result and the equation.
print("The number of such elements is the size of the conjugacy class of a Coxeter element.")
print(f"The calculation is the group order |H3| divided by the Coxeter number h.")
print(f"The order of H3 is the product of its degrees: {degrees[0]} * {degrees[1]} * {degrees[2]} = {group_order}")
print(f"The Coxeter number h for H3 is: {coxeter_number}")
print("The final equation is:")
print(f"{group_order} / {coxeter_number} = {num_elements}")