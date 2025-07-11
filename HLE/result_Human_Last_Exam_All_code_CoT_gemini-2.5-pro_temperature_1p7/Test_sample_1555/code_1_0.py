# In the large-N limit of the CP(N-1) sigma model, the mass spectrum of
# the bound states (excitations) becomes harmonic. The mass of the k-th
# excitation, M_k, is proportional to its quantum number k.

# Define the quantum number for the lightest excitation.
k_lightest = 1

# Define the quantum number for the subsequent higher excitation.
k_subsequent = 2

# The mass ratio is the ratio of their respective quantum numbers.
# Mass Ratio = M_subsequent / M_lightest
# This simplifies to k_subsequent / k_lightest because the proportionality
# constant cancels out.
mass_ratio = k_subsequent / k_lightest

# Print the final equation with all the numbers involved.
print(f"The leading-order asymptotic mass ratio is given by the equation: {k_subsequent} / {k_lightest} = {int(mass_ratio)}")