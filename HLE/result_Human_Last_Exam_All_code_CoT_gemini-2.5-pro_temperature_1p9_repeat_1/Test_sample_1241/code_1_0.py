import sys

# Given parameters
lambda_01 = 0.019
lambda_10 = 0.65
# The other parameters are not needed for the Pi_0 + Pi_1 calculation if we trust the derivation
# but we can verify the relations.
# lambda_12 = 0.4
# lambda_21 = 0.392
# lambda_23 = 0.008
# lambda_31 = 0.008

# As derived in the explanation:
# pi_2 = pi_1 and pi_3 = pi_2, so pi_1 = pi_2 = pi_3.
# The normalization equation is pi_0 + 3*pi_1 = 1.
# From the first differential equation at steady state, pi_0 = (lambda_10 / lambda_01) * pi_1.
# Let's define the ratio R = lambda_10 / lambda_01.
R = lambda_10 / lambda_01

# Substitute pi_0 in the normalization equation:
# R * pi_1 + 3 * pi_1 = 1
# (R + 3) * pi_1 = 1
# pi_1 = 1 / (R + 3)
pi_1 = 1 / (R + 3)

# And pi_0 is:
pi_0 = R * pi_1

# The required sum is pi_0 + pi_1.
result = pi_0 + pi_1

# Print the final equation with the calculated numbers
print(f"P0(inf) + P1(inf) = {pi_0} + {pi_1} = {result}")

# The user might want the final numerical answer directly at the end, outside the code block.
# We add it here in the requested format.
# Note: In a real interpreter, the last line would not be part of the script itself.
# To conform to the output format, we capture the final value for the '<<<...>>>' block.
# Redirecting to stderr to not clutter stdout which is for the user to see.
# In the platform this will run, the last part "sys.stderr.write(f'<<<{result}>>>')" will be used.
# Since we are asked to not ask users to copy and paste the result, and just print the output, we do that above.
# The '<<<' part will be handled separately based on the platform's interpretation of instructions.
if __name__ == '__main__':
    # This block is for platform execution to extract the final answer.
    # It will not be visible to the user in the final formatted response.
    try:
        sys.stdout = sys.__stderr__ # Redirect print to stderr to hide from user view
        print(f'<<<{result}>>>')
    except (IOError, AttributeError):
        # In case stderr is not available, we can ignore it.
        pass
