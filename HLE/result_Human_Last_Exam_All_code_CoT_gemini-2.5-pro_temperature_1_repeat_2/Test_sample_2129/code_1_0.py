import numpy as np
from scipy.special import j0, j1
from scipy.optimize import root_scalar, minimize_scalar
from scipy.integrate import quad

# Step 1: Determine parameters a and lambda

def count_extrema(n):
    """Counts the number of extrema for y2(x) for a given n."""
    # Extrema are roots of y2'(x) = 0, which is J0(2*sqrt(x)) = (20/n)*x
    # Let u = 2*sqrt(x), so x = u^2/4. The equation becomes J0(u) = (5/n)*u^2
    # We need to find the number of positive roots u.
    
    def equation_to_solve(u):
        if u == 0:
            return 1.0 # J0(0)=1, (5/n)*0^2=0 -> not a root for u>0
        return j0(u) - (5.0 / n) * u**2

    u_max = 50 # Search up to u=50 is sufficient
    roots = set()
    
    u_points = np.linspace(0.01, u_max, 5000)
    for i in range(len(u_points) - 1):
        u1, u2 = u_points[i], u_points[i+1]
        # Check if a root is bracketed
        if equation_to_solve(u1) * equation_to_solve(u2) < 0:
            try:
                sol = root_scalar(equation_to_solve, bracket=[u1, u2])
                # round to avoid finding same root multiple times
                roots.add(round(sol.root, 5))
            except (RuntimeError, ValueError):
                continue
    return len(roots)

# For n=10000, we find a
a = count_extrema(10000)

# For n=-2000, we find lambda
lmbda = count_extrema(-2000)

# Step 2: Determine N
# Assumption: y1(x) = exp(-(lmbda-1)x/2) = exp(-x)
# Intersection: exp(-x) = sqrt(x)*J1(2*sqrt(x)) - (10/n)*x^2

# Case n < 0: at most one intersection
# Condition for one intersection is h'(x) <= 0 for all x > 0
# h'(x) = -e^(-x) - J0(2*sqrt(x)) - (20/|n|)x <= 0
# This requires |n| <= 20 / max((-e^(-x) - J0(2*sqrt(x)))/x)
def q_neg(x):
    return (-np.exp(-x) - j0(2 * np.sqrt(x))) / x

res_neg = minimize_scalar(lambda x: -q_neg(x), bounds=(0.1, 10), method='bounded')
M1 = -res_neg.fun
N_neg = int(20 / M1)

# Case n > 0: at most one intersection (0 or 1)
# Condition is min(h(x)) >= 0
# h(x) = exp(-x) - sqrt(x)*J1(2*sqrt(x)) + (10/n)*x^2 >= 0
# This requires n <= 10 / max((-sqrt(x)*J1(2*sqrt(x)) + exp(-x))/x^2)
def q_pos(x):
    if x < 1e-6: # Taylor expansion for x -> 0 to avoid 0/0
        return 1.0 - 2.0*x + (7.0/6.0)*x**2
    return (np.exp(-x) - np.sqrt(x) * j1(2 * np.sqrt(x))) / x**2

res_pos = minimize_scalar(lambda x: -q_pos(x), bounds=(0.1, 20), method='bounded')
M2 = -res_pos.fun
N_pos = int(10 / M2)

N = N_pos + N_neg

# Step 3: Calculate the final value
# (N + lambda) * (y3(x0))^(lambda/a)

x0 = (np.pi / lmbda)**lmbda
C = (a - lmbda) / (lmbda**a)

# y_2s'(x) for n = a*lambda
def y2s_prime(t):
    # n = a*lmbda = 5*3=15. y_2s'(x) = J0(2*sqrt(x)) - (20/15)*x
    return j0(2 * np.sqrt(t)) - (4.0/3.0) * t

# Integrand for y3(x0) using the definition of the fractional integral I^(1/2)
def integrand_y3(t, x_val):
    # This function calculates the integrand for I^(1/2)[y2s_prime]
    # The singularity at t=x_val is handled by the quad integrator
    return y2s_prime(t) / np.sqrt(x_val - t)

# Calculate y3(x0) = -C * I^(1/2)[y2s'(t)](x0)
# I^(1/2)[f] = 1/Gamma(1/2) * integral(...) = 1/sqrt(pi) * integral(...)
integral_val, _ = quad(integrand_y3, 0, x0, args=(x0,))
y3_x0 = -C * (1 / np.sqrt(np.pi)) * integral_val

# Final result
final_result = (N + lmbda) * (y3_x0)**(lmbda / a)

print(f"The parameters are a = {a} and lambda = {lmbda}.")
print(f"The number of integers N for at most one intersection is {N}.")
print(f"The value x0 is ({np.pi:.4f}/{lmbda})^{lmbda} = {x0:.4f}.")
print(f"The coefficient C for the fractional DE is ({a}-{lmbda})/({lmbda}^{a}) = {C:.4f}.")
print(f"The solution y3 at x0 is y3({x0:.4f}) = {y3_x0:.4f}.")
print(f"The final expression is (N + lambda) * (y3(x0))^(lambda/a).")
print(f"({N} + {lmbda}) * ({y3_x0:.4f})^({lmbda}/{a}) = {final_result:.4f}")
print(f"<<<{final_result}>>>")