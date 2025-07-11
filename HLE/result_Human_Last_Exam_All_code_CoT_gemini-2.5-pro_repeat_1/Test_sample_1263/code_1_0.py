# The user wants to know if the calculation is feasible on the Titan architecture.
# Based on the analysis, the architecture's 4-bit constraint on integers (0-15)
# makes it impossible to compute intermediate products required for the calculation.
# For instance, calculating the shell's mass requires multiplying a term for density (proportional to 3)
# by a term for volume (proportional to 8), resulting in an intermediate value of 24, which exceeds the
# maximum representable value of 15. This issue appears in multiple, unavoidable steps of the calculation.
# Thus, the task cannot be completed on Titan.

# The final answer is N0.
print("N0")