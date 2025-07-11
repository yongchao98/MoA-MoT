import sympy

# Let's represent the coordinates and the 1-form symbolically.
x, y = sympy.symbols('x y')
# Let's define a 1-form eta = g(y) dx on a 2-manifold (like a torus or plane)
# According to our argument, if g(y) is not constant, its group of symmetries is not transitive.
# Let's choose g(y) = cos(y), which is not constant.
g = sympy.cos(y)
eta_dx = g
eta_dy = 0 # Our form is eta = g(y) dx

# The exterior derivative d(eta) is (d(eta_dy)/dx - d(eta_dx)/dy) dx wedge dy
# d(eta_dx)/dy is the derivative of g with respect to y
d_eta_dx_dy = sympy.diff(eta_dx, y)

# We check if d(eta) is zero.
is_d_eta_zero = (d_eta_dx_dy == 0)

# Now, let's consider a translation in the y-direction, F(x, y) = (x, y + c)
c = sympy.symbols('c')
y_new = y + c

# The new form F*(eta) would be g(y+c) dx.
# For F to be a symmetry, we need g(y+c) = g(y).
# If g is not constant, this won't hold for all c.
# This indicates the group of symmetries is not transitive.
# Let's check this for our g.
is_g_invariant = sympy.simplify(sympy.cos(y+c) - sympy.cos(y)) == 0


# On the other hand, if g is a constant, say g=k
k = sympy.symbols('k')
eta_const_dx = k
d_eta_const_dx_dy = sympy.diff(eta_const_dx, y)
is_d_eta_const_zero = (d_eta_const_dx_dy == 0)

# And the invariance under y-translation: k = k, which is always true.
# Invariance under x-translation is also trivial.
# So eta = k dx has all translations as symmetries, which form a transitive group.
# And for this form, d(eta) = 0.

# The argument concludes that for the symmetry group to be transitive, g must be constant, which implies d(eta)=0.
# The problem statement tells us the symmetry group IS transitive.
# Therefore, d(eta) must be 0.
final_conclusion = "d(eta) = 0"

print("Step 1: Define a generic 1-form eta = g(x,y) * alpha, where alpha is a closed form defining a homogeneous foliation.")
print("In suitable coordinates, this is eta = g(y) dx.")
print("Step 2: Calculate d(eta) for eta = g(y) dx. This is d(g)/dy dy^dx.")
print("So, d(eta) = 0 if and only if g(y) is a constant.")
print("Step 3: Analyze the symmetry condition. F*(eta) = eta.")
print("For a y-translation F(x,y) = (x, y+c), F*(eta) = g(y+c) dx.")
print("Invariance requires g(y+c) = g(y) for all y and c.")
print("Step 4: The problem states the group of symmetries is transitive. This requires invariance under all translations.")
print("This forces the function g to be a constant.")
print("Step 5: Since g must be a constant, it follows from Step 2 that d(eta) must be 0.")
print("This argument applies to the torus, cylinder, and the plane.")
print("Therefore, it is necessary in any case that d(eta) = 0.")

