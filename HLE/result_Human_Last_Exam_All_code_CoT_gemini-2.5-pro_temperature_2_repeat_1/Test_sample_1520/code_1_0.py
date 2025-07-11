def su_n_generators(n):
  """Calculates the number of generators for an SU(N) group."""
  return n**2 - 1

# Step 1: Calculate the generators for the initial group G = SU(3)
g_group_n = 3
num_g_generators = su_n_generators(g_group_n)
print(f"The initial symmetry group is SU({g_group_n}).")
print(f"Number of generators for SU({g_group_n}) = {g_group_n}^2 - 1 = {num_g_generators}")
print("-" * 30)

# Step 2: Calculate the generators for the residual group H = SU(2) x U(1)
h_group_su_n = 2
h_group_u1_n = 1
num_h_su_generators = su_n_generators(h_group_su_n)
num_h_u1_generators = 1
num_h_generators = num_h_su_generators + num_h_u1_generators
print(f"The residual symmetry group is SU({h_group_su_n}) x U({h_group_u1_n}).")
print(f"Number of unbroken SU({h_group_su_n}) generators = {h_group_su_n}^2 - 1 = {num_h_su_generators}")
print(f"Number of unbroken U({h_group_u1_n}) generators = {num_h_u1_generators}")
print(f"Total number of unbroken generators = {num_h_su_generators} + {num_h_u1_generators} = {num_h_generators}")
print("-" * 30)

# Step 3: Calculate the number of broken generators
num_broken_generators = num_g_generators - num_h_generators
print("The number of broken generators is the difference between the total and unbroken generators.")
print(f"Number of broken generators = (Total generators) - (Unbroken generators)")
print(f"The final equation is: {num_g_generators} - {num_h_generators} = {num_broken_generators}")
print("-" * 30)

# Step 4: Relate to physical consequences
print("In a non-Abelian gauge theory, each broken generator corresponds to a gauge boson that becomes massive via the Higgs mechanism.")
print(f"Therefore, there are {num_broken_generators} massive gauge bosons.")
print("This matches choice E.")
