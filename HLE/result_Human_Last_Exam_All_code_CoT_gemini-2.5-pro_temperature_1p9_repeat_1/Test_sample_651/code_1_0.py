import math

# The problem asks for the limit of M(theta) as theta -> 0.
# Our analysis shows that for any theta > 0, the supremum M(theta) is pi.
# M(theta) is a constant function for theta > 0.
# Let's outline the reasoning for M(theta) = pi:
# 1. The angle 'alpha' is the angle between the trajectory vector (v) and the inner normal vector (n_in).
#    This is calculated as arccos((v . n_in) / (|v|*|n_in|)), which ranges from 0 to pi.
# 2. To maximize this angle, we seek a trajectory where v and n_in are anti-parallel,
#    which would make the angle between them pi.
# 3. The inner normal to side A has an angle of (pi + theta/2).
# 4. The opposite direction thus has an angle of theta/2.
# 5. We showed that a valid trajectory exists with a direction angle of precisely theta/2.
#    This can be constructed by starting at x = exp(i*theta/2) and hitting side A at y = (5/2)*(1 + exp(i*theta)).
#    The vector y-x has exactly the direction theta/2.
# 6. For this trajectory, the angle alpha is pi. Since pi is the maximum possible angle between two vectors,
#    the supremum M(theta) must be pi for any theta > 0.

# The limit of a constant function is the constant itself.
limit_value = math.pi

# The final answer is pi.
# "Remember in the final code you still need to output each number in the final equation!"
# The final result is simply the number pi.
print(limit_value)