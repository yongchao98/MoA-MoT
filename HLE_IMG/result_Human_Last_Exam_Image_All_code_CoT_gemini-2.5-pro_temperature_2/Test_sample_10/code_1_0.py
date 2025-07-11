import numpy as np

# Step 1 & 2: Define the problem based on the number of structural units (n) and known inversion barriers (B).
# We have B(n) = a*n^2 + b*n
# n=1: a*(1**2) + b*1 = 10  =>   a +  b = 10
# n=2: a*(2**2) + b*2 = 49  =>  4a + 2b = 49

# Step 3: Solve the system of linear equations to find a and b.
# We can represent the system as Ax = y, where x = [a, b]
A = np.array([[1, 1], [4, 2]])
y = np.array([10, 49])

# Solve for x = [a, b]
try:
    a, b = np.linalg.solve(A, y)
except np.linalg.LinAlgError:
    # Manual solving if numpy is not available or fails
    # From a + b = 10, we get b = 10 - a
    # Substitute into 4a + 2b = 49:
    # 4a + 2*(10 - a) = 49
    # 4a + 20 - 2a = 49
    # 2a = 29 => a = 14.5
    # b = 10 - 14.5 = -4.5
    a = 14.5
    b = -4.5
    
print("Based on the provided data, we've modeled the inversion barrier B as a function of n units.")
print("The model is a quadratic equation: B(n) = a*n^2 + b*n")
print("\nSolving the system of equations for the coefficients a and b:")
print(f"a = {a}")
print(f"b = {b}")
print(f"So, the predicted relationship is: B(n) = {a}*n^2 + {b}*n")

# Step 4: Predict the barrier for the third molecule, where n=3.
n3 = 3
barrier_n3 = a * (n3**2) + b * n3

print(f"\nNow we calculate the barrier for triacenaphtho[...]triphenylene, which has n = {n3}:")
print(f"B({n3}) = {a} * ({n3}^2) + ({b}) * {n3}")
print(f"B({n3}) = {a} * {n3**2} + ({b}) * {n3}")
print(f"B({n3}) = {a * (n3**2)} - {abs(b * n3)}")
print(f"B({n3}) = {int(round(barrier_n3))}")

print(f"\nThe predicted inversion barrier for triacenaphtho[3,2,1,8-cdefg:3',2',1',8'-ijklm:3'',2'',1'',8''-opqra]triphenylene is {int(round(barrier_n3))} kcal/mol.")
<<<117>>>