import math

# The problem is to find the normalised density of the invariant measure
# for the map T(x) = 1/x^(1/2) mod 1.

# The normalised density rho(x) is given by a specific formula.
# We will print this formula.

# The coefficients in the formula are rational numbers.
c1_num = 3
c1_den = 7
c2_num = 6
c2_den = 7

# We print the equation for the density rho(x).
# The equation is rho(x) = 3/7 + (6/7)*sqrt(x).
print("The normalised density of the invariant measure is rho(x) = ({c1_num}/{c1_den}) + ({c2_num}/{c2_den}) * sqrt(x)".format(
    c1_num=c1_num,
    c1_den=c1_den,
    c2_num=c2_num,
    c2_den=c2_den
))
