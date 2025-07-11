# M: number of spin-degenerate edge states
M = 4
# N: number of reflected spin-degenerate edge states
N = 1

print(f"Calculating the four-terminal conductance G_12,34 for a quantum Hall device.")
print(f"Number of spin-degenerate edge states (M): {M}")
print(f"Number of reflected states (N): {N}\n")

# The derived formula for the conductance G_12,34 is:
# G = 2 * (e^2/h) * (M * (M - N)) / N

print("The general formula for the conductance is: G_12,34 = 2 * (e^2/h) * [M * (M - N)] / N")

if N <= 0 or N > M:
    print("\nError: N must be a positive integer and cannot be greater than M.")
    if N <= 0:
        print("A value of N=0 would imply V_34 = 0, leading to an infinite conductance, which is an ideal limit.")
    if N > M:
        print("The number of reflected states (N) cannot exceed the total number of states (M).")
else:
    # Calculate the numerical coefficient for the conductance quantum e^2/h
    coefficient = 2 * M * (M - N) / N
    
    # Print the result, showing each number in the final equation
    print("\nSubstituting the given values into the formula:")
    print(f"G_12,34 = 2 * (e^2/h) * ({M} * ({M} - {N})) / {N}")
    print(f"G_12,34 = 2 * (e^2/h) * ({M} * {M-N}) / {N}")
    print(f"G_12,34 = 2 * (e^2/h) * ({M * (M-N)}) / {N}")
    
    print("\nThe final result is:")
    print(f"G_12,34 = {coefficient} * e^2/h")
