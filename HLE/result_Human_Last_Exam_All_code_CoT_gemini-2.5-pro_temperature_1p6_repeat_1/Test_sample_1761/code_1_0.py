import sympy

# n is a symbol representing the dimension of the projective space
n = sympy.Symbol('n')

print("To find the complex dimension of the space of global sections of the sheaf Omega^1_{P^n_C} tensor O_{P^n_C}(2),")
print("denoted as h^0(P^n, Omega^1(2)), we follow these steps:")

print("\nStep 1: Use the twisted Euler sequence on P^n.")
print("The standard Euler sequence is: 0 -> Omega^1 -> O(-1)^(n+1) -> O -> 0")
print("Twisting this sequence by the line bundle O(2), we get the short exact sequence:")
print("0 -> Omega^1(2) -> O(1)^(n+1) -> O(2) -> 0")

print("\nStep 2: Consider the associated long exact sequence in cohomology.")
print("This sequence begins with:")
print("0 -> H^0(P^n, Omega^1(2)) -> H^0(P^n, O(1)^(n+1)) -> H^0(P^n, O(2)) -> ...")
print("The map H^0(P^n, O(1)^(n+1)) -> H^0(P^n, O(2)) can be shown to be surjective.")
print("This implies that H^0(P^n, Omega^1(2)) is the kernel of this map.")
print("Therefore, the dimension is given by the equation: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2)).")

# Formulas for cohomology dimensions are expressed using binomial coefficients.
# h^0(P^n, O(k)) is the number of homogeneous polynomials of degree k in n+1 variables,
# which is given by the binomial coefficient C(n+k, k).
h0_O1 = sympy.binomial(n + 1, 1)
h0_O1_oplus = (n + 1) * h0_O1

h0_O2 = sympy.binomial(n + 2, 2)

# Simplify the expressions for printing
h0_O1_oplus_expr = sympy.simplify(h0_O1_oplus)
h0_O2_expr = sympy.simplify(h0_O2)

print("\nStep 3: Substitute the known formulas for the dimensions.")
print("h^0(O(1)^(n+1)) = (n+1) * C(n+1, 1)")
print(f"h^0(O(2))       = C(n+2, 2)")

result_expr = h0_O1_oplus_expr - h0_O2_expr
simplified_result = sympy.simplify(result_expr)

print("\nStep 4: Compute the final dimension by substituting the formulas and simplifying.")
print("The final equation is:")
# The following print statement shows the equation with each "number" (symbolic expression)
print(f"h^0(Omega^1(2)) = ({h0_O1_oplus_expr}) - ({h0_O2_expr})")
print(f"               = {simplified_result}")

print("\nThe final formula for the dimension is:")
print(simplified_result)