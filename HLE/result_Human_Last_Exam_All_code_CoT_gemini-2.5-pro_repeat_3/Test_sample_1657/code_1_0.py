# The problem is unsolvable because the initial distance to Pandora is a required variable for the calculation, but it is not provided.
#
# Let D be the initial distance to Pandora.
# Let v_p(t) be the velocity of Pioneer at time t.
# Let v_r be the recession velocity of Pandora.
#
# The position of Pioneer at time t is: X_p(t) = integral from 0 to t of v_p(t') dt'
# The position of Pandora at time t is: X_pan(t) = D + v_r * t
#
# To find the arrival time T, we must solve for T where X_p(T) = X_pan(T).
# integral from 0 to T of v_p(t') dt' = D + v_r * T
#
# As you can see, the arrival time T is a function of the initial distance D.
# Since D is unknown, T cannot be calculated.
# The problem explicitly states to output "0:0" if a program cannot be written to find the answers.

print("0:0")