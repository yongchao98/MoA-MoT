# Step 1: Define the parameters for the torus knot.
# The knot is identified as the torus knot T(p,q) from the grid diagram.
# The grid size n corresponds to p.
p = 5
# The permutation of the 'x' markers determines q.
# From the problem description, sigma_x = [4, 5, 1, 2, 3], which corresponds to T(5,3).
q = 3

# Step 2: Apply the formula for the maximal Thurston-Bennequin number of T(p,q).
# The formula is tb_max = p*q - p - q.
maximal_tb = p * q - p - q

# Step 3: Print the calculation and the result.
print(f"The grid diagram corresponds to the torus knot T({p},{q}).")
print("The formula for the maximal Thurston-Bennequin number is p*q - p - q.")
print(f"Calculation: {p} * {q} - {p} - {q} = {maximal_tb}")
print(f"The maximal Thurston-Bennequin number is: {maximal_tb}")
