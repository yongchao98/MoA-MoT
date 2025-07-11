def get_su_generators(N):
    """Calculates the number of generators for the SU(N) group."""
    return N**2 - 1

def get_u1_generators():
    """Returns the number of generators for the U(1) group."""
    return 1

# 1. Define the initial and final groups
initial_group_n = 3
final_group_su_n = 2
final_group_has_u1 = True

print(f"The spontaneous symmetry breaking is SU({initial_group_n}) -> SU({final_group_su_n}) x U(1).")
print("-" * 50)

# 2. Calculate generators for the initial group G = SU(3)
dim_G = get_su_generators(initial_group_n)
print(f"The number of generators for the initial group G=SU({initial_group_n}) is:")
print(f"dim(G) = {initial_group_n}^2 - 1 = {dim_G}")
print("-" * 50)

# 3. Calculate generators for the residual group H = SU(2) x U(1)
dim_H_su2 = get_su_generators(final_group_su_n)
dim_H_u1 = get_u1_generators()
dim_H = dim_H_su2 + dim_H_u1

print("The number of generators for the residual (unbroken) group H=SU(2)xU(1) is:")
print(f"dim(SU({final_group_su_n})) = {final_group_su_n}^2 - 1 = {dim_H_su2}")
print(f"dim(U(1)) = {dim_H_u1}")
print(f"Total unbroken generators, dim(H) = {dim_H_su2} + {dim_H_u1} = {dim_H}")
print("-" * 50)

# 4. Calculate the number of broken generators
num_broken_generators = dim_G - dim_H
print("The number of broken generators is the difference dim(G) - dim(H).")
print(f"Number of broken generators = {dim_G} - {dim_H} = {num_broken_generators}")
print("-" * 50)

# 5. Relate broken generators to massive gauge bosons
print("In a non-Abelian gauge theory, the Higgs mechanism makes one gauge boson massive for each broken generator.")
print(f"Therefore, the number of massive gauge bosons is equal to the number of broken generators.")
print(f"Result: There are {num_broken_generators} massive gauge bosons.")
print("-" * 50)
print("This corresponds to option E.")
