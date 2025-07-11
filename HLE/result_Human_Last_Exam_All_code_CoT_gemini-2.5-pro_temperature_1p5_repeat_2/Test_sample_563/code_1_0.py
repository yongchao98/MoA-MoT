# This script provides the solution to a problem whose answer lies in
# established mathematical literature on the theory of Riemann surfaces.
# The numbers of isomorphism classes of automorphism groups for surfaces
# of low genus have been determined by extensive classification work.

# The number of isomorphism classes of automorphism groups for genus g=2.
num_groups_g2 = 13

# The number of isomorphism classes of automorphism groups for genus g=3.
num_groups_g3 = 19

# The number of isomorphism classes of automorphism groups for genus g=4.
num_groups_g4 = 36

# The final answer is a list of these three numbers.
# The prompt requests to "output each number in the final equation".
# We interpret this as showing the individual components of the final answer.
print("The number of isomorphism classes of automorphism groups are as follows:")
print(f"For genus g=2: {num_groups_g2}")
print(f"For genus g=3: {num_groups_g3}")
print(f"For genus g=4: {num_groups_g4}")

# The final result in the requested list format.
result = [num_groups_g2, num_groups_g3, num_groups_g4]
print("\nThe final result as a list for g=[2,3,4] is:")
print(result)