# The equation of the curve is y^2 = f(x), where
# f(x) = 8*x**5 + 4*x**4 + 4*x**3 + 1*x**2 + 8*x
#
# Step 1: Find the genus of the curve.
# The degree 'd' of the polynomial f(x) is 5.
d = 5

# The genus 'g' of a hyperelliptic curve y^2 = f(x) is given by g = floor((d-1)/2).
# Since d=5 is odd, g = (d-1)/2.
g = (d - 1) / 2

# Step 2: Relate the genus to the stable reduction.
# The arithmetic genus of a curve is preserved under stable reduction.
# Let's assume the stable reduction of this curve is a single irreducible rational curve.
# A rational curve has geometric genus 0.
# The arithmetic genus of an irreducible curve is the sum of its geometric genus
# and its number of double points (delta).
# p_a_reduced = p_g_reduced + delta
#
# Step 3: Calculate the number of double points.
# By equating the genus of the original curve with the arithmetic genus of its stable reduction:
# g = p_a_reduced = 0 + delta
# Thus, delta = g.
num_double_points = g

# The final equation determines the number of double points based on the degree of the polynomial.
print("The number of double points is found by first calculating the genus of the curve.")
print("The calculation for the number of double points is:")
print(f"{int(num_double_points)} = ({d} - 1) / 2")
