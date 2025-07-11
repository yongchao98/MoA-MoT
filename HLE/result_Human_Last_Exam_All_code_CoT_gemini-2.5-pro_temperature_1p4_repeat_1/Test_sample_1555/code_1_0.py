import math

# Step 1: Explain the physical context and the mass formula.
print("In the two-dimensional CP(N-1) model, the mass spectrum of the tower of solitonic excitations (or kinks) can be determined using semi-classical methods.")
print("The mass of the k-th excitation, M_k, for k=1, 2, ..., N-1, is given by the formula:")
print("  M_k = C * sin(k * pi / N)")
print("where C is a constant proportional to the dynamically generated mass scale of the theory.")
print("-" * 50)

# Step 2: Identify the specific masses for the ratio.
print("The lightest solitonic excitation corresponds to k=1:")
print("  M_1 = C * sin(1 * pi / N)")
print("\nThe subsequent higher excitation corresponds to k=2:")
print("  M_2 = C * sin(2 * pi / N)")
print("-" * 50)

# Step 3: Formulate and simplify the mass ratio.
print("The ratio of these two masses is M_2 / M_1:")
print("  Ratio = (C * sin(2 * pi / N)) / (C * sin(pi / N))")
print("  Ratio = sin(2 * pi / N) / sin(pi / N)")
print("\nUsing the trigonometric double-angle identity sin(2x) = 2*sin(x)*cos(x), we simplify the ratio:")
print("  Ratio = (2 * sin(pi / N) * cos(pi / N)) / sin(pi / N)")
print("  Ratio = 2 * cos(pi / N)")
print("-" * 50)

# Step 4: Evaluate the asymptotic limit as N -> infinity.
print("The problem asks for the leading-order asymptotic mass ratio as N approaches infinity.")
print("In this limit (N -> infinity), the argument of the cosine term (pi / N) approaches 0.")
print("Therefore, we evaluate the limit of cos(pi / N):")
print("  lim (N->inf) cos(pi / N) = cos(0)")
print("-" * 50)

# Step 5: Calculate the final numerical answer.
print("The final calculation for the ratio is:")
factor = 2
cos_of_zero = math.cos(0)
final_ratio = factor * cos_of_zero

print(f"  {factor} * cos(0) = {factor} * {cos_of_zero} = {final_ratio}")
print("\nThus, the leading-order asymptotic mass ratio is 2.")
<<<2>>>