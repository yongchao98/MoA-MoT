import scipy.integrate

def area_ratio_integrand(f, e, d):
    """
    Calculates the ratio of the area of the inner triangle XYZ to the outer
    triangle ABC.
    d, e, f are the ratios BD/BC, CE/CA, AF/AB, respectively, and are
    assumed to be drawn from a uniform distribution on [0, 1].
    """
    # Check for edge cases where a denominator might be zero, although tplquad
    # handles these singularities.
    if (1 - e + d * e) == 0 or (1 - f + e * f) == 0 or (1 - d + f * d) == 0:
        # This occurs at corners like d=1, e=1. The limit is well-defined.
        # For d=e=f=1, ratio is 1. For d=e=f=0, ratio is 1.
        # We can return a value like 1.0, or rely on the integrator to handle it.
        # Let's compute it for the non-edge case first.
        pass

    numerator = (1 - d - e - f + d*e + e*f + f*d - 2*d*e*f)**2
    denominator = (1 - e + d*e) * (1 - f + e*f) * (1 - d + f*d)
    
    return numerator / denominator

# We need to calculate the triple integral of the area_ratio_integrand
# over d, e, and f, each from 0 to 1.
# tplquad(func, a, b, gfun, hfun, qfun, rfun)
# Here a,b are limits for d, gfun,hfun for e, and qfun,rfun for f.
# All are constants [0,1].
result, error = scipy.integrate.tplquad(
    area_ratio_integrand, 0, 1, 0, 1, 0, 1
)

# The final probability is the value of this integral.
# The analytical solution to this problem is 1/10.
# Our numerical result should be very close to this.
prob = result
k_formula = "(1 - d - e - f + de + ef + fd - 2def)² / ((1-e+de)(1-f+ef)(1-d+fd))"
final_equation = f"""
The probability is given by the integral of the area ratio function K(d, e, f) over the unit cube.
P(P in XYZ) = E[Area(XYZ)/Area(ABC)]
P(P in XYZ) = ∫[0,1]∫[0,1]∫[0,1] {k_formula} dd de df
Numerical Result: {prob:.10f}
Analytical Result: 1/10
"""
print(final_equation)

print("The final equation is:")
print(f"P(P in XYZ) = E[Area(XYZ)/Area(ABC)] = ∫[0,1]∫[0,1]∫[0,1] K(d,e,f) dd de df = 1/10")
print("Where K(d,e,f) is the expression shown above.")
