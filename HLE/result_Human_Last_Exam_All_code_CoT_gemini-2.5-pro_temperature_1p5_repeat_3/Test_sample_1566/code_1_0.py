# The problem asks for the number of topologically distinct continua satisfying two properties.
# Through topological analysis, we deduce the following:

# 1. The properties imply the continuum must be a "chainable continuum"
#    with exactly two endpoints.

# 2. The properties also imply that the auto-homeomorphism group of the continuum
#    must be transitive on its interior points.

# 3. A key theorem in continuum theory states that the only such continuum
#    is the arc (a space homeomorphic to [0,1]). The arc is a Peano continuum.
num_peano_continua = 1

# 4. The theorem proves the uniqueness of the arc, meaning no non-Peano continua
#    (which are not locally connected) can satisfy the conditions.
num_non_peano_continua = 0

# 5. The total number of distinct continua is the sum of these cases.
total_continua = num_peano_continua + num_non_peano_continua

# The final equation showing how the total is calculated:
print(f"The number of topologically distinct continua is the sum of possibilities:")
print(f"Peano (arc-like): {num_peano_continua}")
print(f"Non-Peano: {num_non_peano_continua}")
print(f"Total: {num_peano_continua} + {num_non_peano_continua} = {total_continua}")
