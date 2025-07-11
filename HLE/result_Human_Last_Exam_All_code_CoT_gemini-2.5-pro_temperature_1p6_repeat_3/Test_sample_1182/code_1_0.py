# The genus of the original curve y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5
# is g = floor((deg(polynomial) - 1) / 2).
# The degree of the polynomial is 5.
g = (5 - 1) // 2

# Through a sequence of transformations, we find a regular model of the curve
# whose special fiber modulo 2 is a smooth curve (a conic).
# A smooth conic has genus 0.
gs = 0

# The number of double points (delta) in the stable reduction is the
# difference between the genus of the original curve (g) and the genus
# of the special fiber of its regular model (gs).
delta = g - gs

print(f"The genus of the original curve is g = {g}.")
print(f"The genus of the special fiber of the regular model is gs = {gs}.")
print("The number of double points in the stable reduction is given by the formula: delta = g - gs.")
print(f"So, the number of double points is {g} - {gs} = {delta}.")