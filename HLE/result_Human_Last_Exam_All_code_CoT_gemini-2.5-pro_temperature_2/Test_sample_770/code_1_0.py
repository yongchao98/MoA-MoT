import math

print("This script calculates the rank of H^2_c(Y, Q) for the given space Y.")
print("The rank is determined by the number of 'age 1' conjugacy classes of the icosahedral group G = A_5 in SL(3,C).")
print("\nWe check the age for each of the 4 non-trivial conjugacy classes.")
print("The age of a group element g is the sum of the fractional parts of its eigenvalue exponents (the alphas).")
print("For any rotation in SO(3), one eigenvalue is 1, so one alpha is always 0.")
print("-" * 80)

# The non-trivial conjugacy classes of A_5 correspond to rotations.
# The eigenvalues for a rotation by theta are 1, exp(i*theta), exp(-i*theta).
# In our convention (0 <= alpha < 1), these are exp(2*pi*i*alpha) with alphas 0, theta/(2*pi), (2*pi-theta)/(2*pi).

# Class 1: Rotations by 2*pi/3 (order 3 elements)
class1_name = "3-cycles (rotation by 2pi/3)"
alpha1_1, alpha1_2, alpha1_3 = 0, 1/3, 2/3
age1 = alpha1_1 + alpha1_2 + alpha1_3
print(f"For the class '{class1_name}':")
print(f"  Eigenvalue exponents (alphas) are {alpha1_1}, {alpha1_2:.3f}, {alpha1_3:.3f}")
print(f"  Age = {alpha1_1} + {alpha1_2:.3f} + {alpha1_3:.3f} = {age1:.1f}. This is an age 1 class.")
print("-" * 80)

# Class 2: Rotations by pi (order 2 elements, double transpositions)
class2_name = "Double transpositions (rotation by pi)"
alpha2_1, alpha2_2, alpha2_3 = 0, 1/2, 1/2
age2 = alpha2_1 + alpha2_2 + alpha2_3
print(f"For the class '{class2_name}':")
print(f"  Eigenvalue exponents (alphas) are {alpha2_1}, {alpha2_2:.3f}, {alpha2_3:.3f}")
print(f"  Age = {alpha2_1} + {alpha2_2:.3f} + {alpha2_3:.3f} = {age2:.1f}. This is an age 1 class.")
print("-" * 80)

# Class 3: Rotations by 2*pi/5 (order 5 elements, one class of 5-cycles)
class3_name = "5-cycles, type 1 (rotation by 2pi/5)"
alpha3_1, alpha3_2, alpha3_3 = 0, 1/5, 4/5
age3 = alpha3_1 + alpha3_2 + alpha3_3
print(f"For the class '{class3_name}':")
print(f"  Eigenvalue exponents (alphas) are {alpha3_1}, {alpha3_2:.3f}, {alpha3_3:.3f}")
print(f"  Age = {alpha3_1} + {alpha3_2:.3f} + {alpha3_3:.3f} = {age3:.1f}. This is an age 1 class.")
print("-" * 80)

# Class 4: Rotations by 4*pi/5 (order 5 elements, the other class of 5-cycles)
class4_name = "5-cycles, type 2 (rotation by 4pi/5)"
alpha4_1, alpha4_2, alpha4_3 = 0, 2/5, 3/5
age4 = alpha4_1 + alpha4_2 + alpha4_3
print(f"For the class '{class4_name}':")
print(f"  Eigenvalue exponents (alphas) are {alpha4_1}, {alpha4_2:.3f}, {alpha4_3:.3f}")
print(f"  Age = {alpha4_1} + {alpha4_2:.3f} + {alpha4_3:.3f} = {age4:.1f}. This is an age 1 class.")
print("-" * 80)

# The identity element has age 0 and is not counted.
count_c1 = 1 if math.isclose(age1, 1) else 0
count_c2 = 1 if math.isclose(age2, 1) else 0
count_c3 = 1 if math.isclose(age3, 1) else 0
count_c4 = 1 if math.isclose(age4, 1) else 0

total_count = count_c1 + count_c2 + count_c3 + count_c4

print("\nFinal Result:")
print("The rank of H^2_c(Y, Q) is the total number of non-trivial conjugacy classes with age 1.")
print("The final equation is the sum of indicators (1 if age=1, 0 otherwise) for these classes:")
print(f"{count_c1} + {count_c2} + {count_c3} + {count_c4} = {total_count}")