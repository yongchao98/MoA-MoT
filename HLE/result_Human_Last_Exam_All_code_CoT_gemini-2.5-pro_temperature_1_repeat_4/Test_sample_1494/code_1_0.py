import sympy

# n_power represents the power of k in the Taylor expansion of the initial data's Fourier transform at k=0.
# For the best-case scenario (zero horizontal mean), the transform starts with k^1, so its square |a_k|^2 starts with k^2.
# We are calculating the norm of the x-derivative, so we have an extra k^2 factor.
# Thus, the integrand for the squared norm is k^2 * |a_k|^2 ~ k^4.
# The general form of the integral is integral(k^(2*n) * exp(-a*k^2)) dk
# For the best-case decay of ||d_x rho||^2, the power is 4, so n=2.

n = 2 # Set n=2 for the k^4 case (zero-mean initial data)

# Define symbols
k = sympy.Symbol('k', positive=True)
a = sympy.Symbol('a', positive=True)

# The integral that determines the decay rate
# We are interested in how this integral scales with 'a'
integral_expr = k**(2*n) * sympy.exp(-a * k**2)

# Compute the definite integral from 0 to oo
# This integral gives the time-dependent coefficient of the decay
result = sympy.integrate(integral_expr, (k, 0, sympy.oo))

print(f"The integral of k^({2*n}) * exp(-a*k^2) from 0 to oo is:")
print(result)
print("-" * 20)
# The parameter 'a' scales linearly with time 't'.
# We need to find the power of 'a' in the result.
# The result is of the form C * a^exponent, where C is a constant.
# We extract the exponent of 'a'.
exponent_of_a = result.as_powers_dict()[a]

print(f"The scaling of the integral with respect to 'a' is: a^({exponent_of_a})")

# The decay of ||d_x rho||^2 is proportional to this scaling. Since a ~ t, ||d_x rho||^2 ~ t^exponent.
# The decay of ||d_x rho|| is the square root of this.
exponent_of_t = exponent_of_a / 2

print(f"Since a is proportional to t, the decay of ||d_x rho||^2 is like t^({exponent_of_a}).")
print(f"Therefore, the best time-decay we can expect for ||d_x rho(t)||_L2 is of the form t^({exponent_of_t}).")

final_decay_power_num, final_decay_power_den = exponent_of_t.as_numer_denom()

print("\nThe final decay is expressed as t^(-5/4).")
print("The exponent is:")
print(f"{final_decay_power_num}/{final_decay_power_den}")
