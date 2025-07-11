# This code must be run in a SageMath environment.

# 1. Define the polynomial ring over the rational numbers.
P.<x> = QQ['x']

# 2. Define the polynomial from the curve's equation.
# y^2 = x^6 + 4x^5 + 6x^4 + 2x^3 + x^2 + 2x + 1
f = x^6 + 4*x^5 + 6*x^4 + 2*x^3 + x^2 + 2*x + 1

# 3. Define the hyperelliptic curve.
# As established in the plan, this is a genus 2 curve.
C = HyperellipticCurve(f)

# 4. Get the Jacobian of the curve.
J = C.jacobian()

# 5. Calculate the conductor of the Jacobian.
conductor = J.conductor()

# 6. Print the conductor.
print(conductor)