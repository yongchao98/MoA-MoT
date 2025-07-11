import math

# Based on the mathematical derivation, for any small theta > 0,
# it's possible to have a trajectory with angle psi = theta/2.
# For such a trajectory, the variable delta = psi - theta/2 becomes 0.
# The angle alpha is given by alpha = arccos(-cos(delta)).

# The supremum of alpha, M(theta), corresponds to the maximum possible value of cos(delta).
# Since delta=0 is achievable, the maximum of cos(delta) is cos(0).
max_cos_delta = math.cos(0)

# The value inside the arccos function is -cos(delta), so its minimum is -1.
value_in_arccos = -max_cos_delta

# The supremum M(theta) is the arccos of this minimum value.
M_theta = math.acos(value_in_arccos)

# The value of M(theta) is constant for all small theta > 0.
# Therefore, its limit as theta -> 0 is this constant value.

print(f"The maximum value for cos(delta) is found to be: {max_cos_delta}")
print(f"This leads to a value of {value_in_arccos} inside the arccos function for alpha.")
print(f"The final equation for the supremum M(theta) is arccos({value_in_arccos}).")
print(f"The result is {M_theta}.")

final_answer = M_theta
# We can represent the final answer using the symbolic constant pi.
print(f"The limit is pi, which is approximately {final_answer}.")
print("<<<pi>>>")