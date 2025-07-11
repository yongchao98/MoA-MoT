# The dimension of the parameter space (number of coefficients a_i).
N = 9

# The minimal decay exponent (slowest decay rate) for the oscillatory integral.
# This corresponds to the case where the phase polynomial is most degenerate,
# effectively reducing to a cubic term in one variable, like x^3.
# The decay rate for an integral with phase lambda * x^k is lambda^(-1/k).
# For k=3, the exponent is 1/3.
sigma_min = 1/3

# The integrability condition for I(a) in L^p(R^N) is p * sigma_min > N.
# The function is not in L^p(R^N) if p * sigma_min <= N.
# We are looking for the largest p that satisfies this condition.
# This occurs at the boundary: p = N / sigma_min.
p = N / sigma_min

print(f"The dimension of the parameter space is N = {N}")
print(f"The minimal decay exponent is sigma_min = 1/{int(1/sigma_min)}")
print(f"The critical value p is given by the equation p = N / sigma_min")
print(f"p = {N} / (1/3) = {int(p)}")
print(f"So, the largest p such that the function I is not in L^p(R^9) is {int(p)}.")

# The final answer in the requested format
# print(f'<<<{int(p)}>>>')