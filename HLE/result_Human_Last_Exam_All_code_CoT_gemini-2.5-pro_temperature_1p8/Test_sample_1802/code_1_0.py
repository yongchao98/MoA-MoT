# The dimension of the qudit system.
# This can be any integer greater than or equal to 2.
# For example, for a qubit d=2, for a qutrit d=3.
d = 3

# Step 1: A Pauli channel is defined by d*d generalized Pauli operators.
num_pauli_operators = d**2

# Step 2: The rank of the channel is the number of linearly independent Kraus operators.
# This corresponds to the number of non-zero probabilities p_jk.
# To maximize the rank, we choose all d*d probabilities to be non-zero.
# For example, p_jk = 1/(d*d) for all j,k. This gives d*d linearly independent Kraus operators.
# Therefore, the maximal rank of the Pauli channel is d*d.
max_rank_Lambda = num_pauli_operators

# Step 3: The rank of a channel is equal to the rank of its complementary channel.
max_rank_complementary = max_rank_Lambda

# Step 4: Print the result and the equation.
print(f"For a qudit of dimension d = {d}, a Pauli channel is characterized by d^2 = {d*d} generalized Pauli operators.")
print("The rank of the channel corresponds to the number of these operators with non-zero probability.")
print("To maximize the rank, we can assign a non-zero probability to all operators.")
print(f"Thus, the maximal rank of the Pauli channel, and its complementary channel, is d^2.")
print(f"Final equation: {d} ^ 2 = {max_rank_complementary}")
