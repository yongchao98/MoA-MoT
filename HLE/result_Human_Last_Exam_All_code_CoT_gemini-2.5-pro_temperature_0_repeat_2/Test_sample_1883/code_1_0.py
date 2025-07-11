# 1. Define the number of pi electrons in the conjugated system.
pi_electrons = 8

# 2. Determine the principal quantum number 'n' of the HOMO for a thermal reaction.
# The HOMO is the orbital with quantum number n = pi_electrons / 2.
n_homo = pi_electrons / 2

# 3. Calculate the number of nodes in the HOMO.
# The number of nodes in a molecular orbital psi_n is n - 1.
num_nodes = n_homo - 1

# 4. The major product (B) is formed via the allowed conrotatory pathway.
# The minor product (A) is formed via a disfavored pathway.
# A simple model predicts the ratio of major to minor product to be N^2,
# where N is the number of nodes in the HOMO.
major_to_minor_ratio_value = num_nodes**2

# 5. Print the final predicted ratio of A to B.
# The equation is Ratio A : B = 1 : N^2
print("The reacting system is an octatetraene with 8 pi electrons.")
print(f"The HOMO for a thermal reaction is psi_{int(n_homo)}.")
print(f"The number of nodes (N) in the HOMO is {int(n_homo)} - 1 = {int(num_nodes)}.")
print("The ratio of the minor product (A) to the major product (B) is predicted by the model 1 : N^2.")
print(f"Therefore, the predicted ratio is:")
print(f"A : B = 1 : {int(num_nodes)}^2")
print(f"A : B = 1 : {int(major_to_minor_ratio_value)}")