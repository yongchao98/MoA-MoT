import numpy

print("We have derived two possible cases for the relationship between a1 and d from the problem conditions.")
print("The script will now analyze each case to find the value of d.")
print("\n------------------------------------")

print("\nStep 1: Analyzing Case 1 where a1 = d.")
print("The condition S_99 - T_99 = 99 leads to the quadratic equation for d:")
# The equation is 50d^2 - d - 51 = 0.
coeffs1 = [50, -1, -51]
print(f"Equation: {coeffs1[0]} * d^2 + ({coeffs1[1]}) * d + ({coeffs1[2]}) = 0")

roots1 = numpy.roots(coeffs1)
print(f"The solutions for d are: {roots1[0]:.4f} and {roots1[1]:.4f}")

# Check for a valid solution (d > 1)
valid_d = None
for r in roots1:
    if r > 1:
        valid_d = r
        print(f"Found a valid solution d = {r:.4f}, which is > 1.")

print("\n------------------------------------")

print("\nStep 2: Analyzing Case 2 where a1 = 2d.")
print("The condition S_99 - T_99 = 99 leads to the quadratic equation for d:")
# The equation is 51d^2 - d - 50 = 0.
coeffs2 = [51, -1, -50]
print(f"Equation: {coeffs2[0]} * d^2 + ({coeffs2[1]}) * d + ({coeffs2[2]}) = 0")
roots2 = numpy.roots(coeffs2)
print(f"The solutions for d are: {roots2[0]:.4f} and {roots2[1]:.4f}")

# Check for a valid solution (d > 1)
is_valid_case2 = False
for r in roots2:
    if r > 1:
        is_valid_case2 = True
if not is_valid_case2:
    print("Neither solution in this case is greater than 1.")
    
print("\n------------------------------------")

# Final Conclusion
print("\nStep 3: Final Conclusion.")
if valid_d is not None:
    print("The only solution satisfying the condition d > 1 comes from Case 1.")
    print(f"The final value of d is {valid_d}")
else:
    print("No valid solution for d was found.")