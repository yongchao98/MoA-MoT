# Set the value for q
q = 997

# --- Part 1: Involutions from g^2 = I ---
# This corresponds to a single class of involutions in PSU(4,q).
# The formula for the number of these involutions is q^4 * (q^2 + 1) * (q^2 - q + 1)
n1 = q**4 * (q**2 + 1) * (q**2 - q + 1)

# --- Part 2: Involutions from g^2 = -I ---
# These arise from elements g in SU(4,q) such that g^2 = -I.
# These elements fall into two conjugacy classes in SU(4,q).

# Number of elements in the first class
n2_class1 = q**2 * (q**3 + 1)

# Number of elements in the second class
n2_class2 = q**4 * (q + 1) * (q**2 - q + 1) * (q**2 + 1)

# The number of resulting involutions in PSU(4,q) is half the total number of these elements.
n2 = (n2_class1 + n2_class2) // 2

# --- Total Number of Involutions ---
total_involutions = n1 + n2

# --- Print the results step-by-step ---
print(f"For q = {q}:")
print("-" * 30)
print("1. Number of involutions from elements g with g^2 = I:")
print(f"N1 = {q}^4 * ({q}^2 + 1) * ({q}^2 - {q} + 1)")
print(f"N1 = {n1}")
print("-" * 30)
print("2. Number of involutions from elements g with g^2 = -I:")
print("This requires counting the elements g and dividing by 2.")
print(f"   Number of elements in class 1: {q}^2 * ({q}^3 + 1) = {n2_class1}")
print(f"   Number of elements in class 2: {q}^4 * ({q}+1) * ({q}^2-{q}+1) * ({q}^2+1) = {n2_class2}")
print(f"N2 = ({n2_class1} + {n2_class2}) / 2")
print(f"N2 = {n2}")
print("-" * 30)
print("Total number of involutions in PSU(4, 997):")
print(f"Total = N1 + N2")
print(f"Total = {n1} + {n2}")
print(f"Total = {total_involutions}")
