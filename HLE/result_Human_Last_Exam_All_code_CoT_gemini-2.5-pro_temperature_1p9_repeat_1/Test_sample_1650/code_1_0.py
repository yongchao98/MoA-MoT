import sympy

# Step 1: Define all necessary mathematical symbols.
# We use 'lam' for lambda to avoid conflict with Python's keyword.
# R is the internuclear distance, and zeta is the effective nuclear charge.
zeta, R, lam, mu, phi = sympy.symbols('zeta R lambda mu phi', real=True, positive=True)

# Step 2: Define the components of the integral.
# The normalized 2s wavefunction is Ψ_2s = N * (2 - ζ*r) * exp(-ζ*r/2).
# The square of the normalization constant is N^2 = ζ^3 / (32*π).
N_sq = zeta**3 / (32 * sympy.pi)

# In elliptical coordinates, the distances from the two nuclei (r_a, r_b) are:
r_a = R/2 * (lam + mu)
r_b = R/2 * (lam - mu)

# The product of the two wavefunctions, Ψ_a * Ψ_b, can be constructed.
# The polynomial part is (2 - ζ*r_a) * (2 - ζ*r_b).
poly_part = (2 - zeta * r_a) * (2 - zeta * r_b)

# The exponential part is exp(-ζ*r_a/2) * exp(-ζ*r_b/2) = exp(-ζ*(r_a+r_b)/2).
# In elliptical coordinates, this simplifies since r_a + r_b = R*λ.
exp_part = sympy.exp(-zeta * R * lam / 2)

# The volume element dτ in elliptical coordinates is (R^3/8) * (λ^2 - μ^2) dλ dμ dφ.
vol_element = (R**3 / 8) * (lam**2 - mu**2)

# Step 3: Combine all parts into the full integrand.
integrand = N_sq * poly_part * exp_part * vol_element

# Step 4: Perform the triple integration symbolically.
# The limits are: φ from 0 to 2π, μ from -1 to 1, λ from 1 to infinity.

# Integration with respect to φ is trivial as the integrand is independent of it.
integral_phi = sympy.integrate(integrand, (phi, 0, 2*sympy.pi))

# Integration with respect to μ.
integral_mu = sympy.integrate(sympy.simplify(integral_phi), (mu, -1, 1))

# Integration with respect to λ.
overlap_integral = sympy.integrate(sympy.simplify(integral_mu), (lam, 1, sympy.oo))

# Step 5: Simplify the final expression and make it more readable.
# The result is typically expressed in terms of a single parameter p = ζ*R.
p = sympy.Symbol('p', positive=True)
final_expression = sympy.simplify(overlap_integral).subs(zeta * R, p)

# Step 6: Print the final result in a clear, formatted equation.
# The result is factored into an exponential and a polynomial part.
final_poly = sympy.collect(sympy.expand(final_expression / sympy.exp(-p/2)), p)

print("The analytical expression for the overlap integral S, with p = ζR, is:")
# The f-string automatically converts the symbolic expression to a string.
# The polynomial part includes all the required numbers from the calculation.
print(f"S(p) = exp(-p/2) * ({final_poly})")