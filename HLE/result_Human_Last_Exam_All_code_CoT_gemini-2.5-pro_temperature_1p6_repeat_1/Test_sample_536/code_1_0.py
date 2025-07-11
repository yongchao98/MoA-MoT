import math

# Given inner products
c_hb = 0.9375  # lim <h_p, b_p>
c_hz = 0.9     # lim <h_p, z_p>

# The problem reduces to solving a quadratic equation for c = lim <b_p, z_p>
# of the form a*c^2 + b*c + c_ = 0

a = 1.0
b = -2 * c_hb * c_hz
c_ = c_hb**2 + c_hz**2 - 1

print("The relationship between the inner products forms a quadratic equation a*c^2 + b*c + c_ = 0 for c = lim <b_p, z_p>.")
print(f"The coefficient a is: {a}")
print(f"The coefficient b is: {b}")
print(f"The coefficient c_ is: {c_}")

# Calculate the discriminant
discriminant = b**2 - 4 * a * c_

if discriminant < 0:
    print("The equation has no real solutions.")
else:
    # Calculate the two possible solutions for c
    sqrt_discriminant = math.sqrt(discriminant)
    sol1 = (-b + sqrt_discriminant) / (2 * a)
    sol2 = (-b - sqrt_discriminant) / (2 * a)

    print(f"\nThe two potential solutions are {sol1} and {sol2}.")
    
    # In this type of problem from random matrix theory, the eigenvector of the
    # sample covariance matrix (h) is typically found to lie geometrically 
    # between the eigenvector of the population covariance matrix (b) and other
    # influential vectors (z). This configuration corresponds to the smaller
    # of the two roots.
    final_answer = min(sol1, sol2)
    print(f"\nThe physically plausible solution is the smaller root.")
    print(f"The final result for lim <b_p, z_p> is: {final_answer}")
    print(f"<<<{final_answer}>>>")
