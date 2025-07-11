import scipy.integrate as integrate
import numpy as np

# The problem is to compute the integral of f(x,y,z) = z^2*(x^2+y^2)
# over a cone with base radius R=3 at z=0 and height H=2.
#
# In cylindrical coordinates:
# f(r,θ,z) = z^2 * r^2
# The volume element dV = r*dr*dθ*dz.
# The integrand becomes (z^2*r^2) * r = z^2*r^3.
#
# The integration limits are:
# - θ from 0 to 2π
# - z from 0 to 2
# - r from 0 to 1.5*(2-z)

# The integral is ∫(θ=0..2π) ∫(z=0..2) ∫(r=0..1.5*(2-z)) z^2*r^3 dr dz dθ.
#
# We use SciPy's tplquad. Its integration order is fixed, so we map our
# variables: outer(x)=θ, middle(y)=z, inner(z)=r.
# The integrand function arguments are in the order (inner, middle, outer),
# so we define it as f(r, z, θ).
integrand = lambda r, z, theta: z**2 * r**3

# Perform the triple integral.
# tplquad(func, outer_min, outer_max, middle_min, middle_max, inner_min, inner_max)
result, error = integrate.tplquad(integrand,
                                  0, 2 * np.pi,             # theta limits
                                  lambda theta: 0,          # z lower limit
                                  lambda theta: 2,          # z upper limit
                                  lambda theta, z: 0,       # r lower limit
                                  lambda theta, z: 1.5 * (2 - z)) # r upper limit

print("The integral of f(x,y)=z^2(x^2+y^2) over the cone is calculated by the following final equation:")
# We print each number in the equation as requested.
print("Integral = ∫(θ from 0 to 2π) dθ ∫(z from 0 to 2) dz ∫(r from 0 to 1.5*(2-z)) (z^2 * r^3) dr")

# The exact analytical result is 108*π/35. We compare our numerical result to it.
exact_value = 108 * np.pi / 35
print(f"\nNumerical Result: {result}")
print(f"Exact Result (108*π/35): {exact_value}")
