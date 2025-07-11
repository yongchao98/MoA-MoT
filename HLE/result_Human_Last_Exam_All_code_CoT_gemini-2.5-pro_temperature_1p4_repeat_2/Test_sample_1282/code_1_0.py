import math

# This script calculates the estimated blow-up time for the given PDE,
# based on the analysis of an ordinary differential equation (ODE) derived from it.

# The enstrophy of the solution, Y(t) = ||∇u(t)||^2_{L^2}, is shown to satisfy
# the differential inequality: dY/dt <= K * Y(t)^3 / (1 + t)^3
# where K is a constant.

# Analysis of the corresponding ODE shows that a blow-up can occur at a finite time T
# if the initial condition Y_0 = Y(0) is sufficiently large, specifically if K * Y_0^2 > 1.
# The blow-up time T is given by the formula:
# T = sqrt( (K * Y_0^2) / (K * Y_0^2 - 1) ) - 1

# Let's assign some example values for the constant K and the initial enstrophy Y_0
# that satisfy the blow-up condition.
K = 1.0
Y_0 = 2.0  # Note: Y_0 is ||∇u_0||^2, not the norm itself.

# The blow-up condition is K * Y_0^2 > 1
blow_up_condition_value = K * Y_0**2

if blow_up_condition_value > 1:
    print(f"For K = {K} and Y_0 = {Y_0}, the blow-up condition is met:")
    print(f"K * Y_0^2 = {K} * {Y_0}^2 = {blow_up_condition_value:.2f}, which is > 1.\n")

    # Calculate the terms in the formula
    numerator = K * Y_0**2
    denominator = K * Y_0**2 - 1
    
    # Calculate the estimated blow-up time
    T_blowup = math.sqrt(numerator / denominator) - 1

    print("The estimated blow-up time T is calculated as follows:")
    print(f"T = sqrt( (K * Y_0^2) / (K * Y_0^2 - 1) ) - 1")
    print(f"T = sqrt( ({K} * {Y_0**2}) / ({K} * {Y_0**2} - 1) ) - 1")
    print(f"T = sqrt( {numerator:.2f} / {denominator:.2f} ) - 1")
    print(f"T = {T_blowup:.4f}")
else:
    print(f"For K = {K} and Y_0 = {Y_0}, the blow-up condition is NOT met:")
    print(f"K * Y_0^2 = {K} * {Y_0}^2 = {blow_up_condition_value:.2f}, which is not > 1.")
    print("A blow-up is not predicted by this model for these values.")
