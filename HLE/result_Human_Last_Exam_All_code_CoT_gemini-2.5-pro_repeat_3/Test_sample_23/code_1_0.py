# Let a and b be positive integers.
# You can change these values to test with different dimensions.
a = 4
b = 5

# The problem asks for the number of "non-admissible" integers k
# in the range 0, 1, ..., a*b. An integer k is admissible if there exists
# an orthogonal basis for the space of a x b matrices where exactly k
# of the basis matrices have rank 1.

# The solution depends on whether the matrix dimensions a and b are greater than 1.
# We analyze two distinct cases.

# Case 1: One of the dimensions is 1 (i.e., the matrices are vectors).
if a == 1 or b == 1:
    # In this case, any non-zero a x b matrix has rank 1.
    # An orthogonal basis for the space M_{a,b}(C) must consist of a*b non-zero matrices.
    # Therefore, all a*b matrices in the basis must have rank 1.
    # This means k = a*b is the only admissible value.
    # The non-admissible integers are all values in {0, 1, ..., a*b} except for a*b.
    # The count of non-admissible integers is thus (a*b + 1) - 1 = a*b.
    result = a * b
    print(f"For a={a}, b={b}, since one dimension is 1, the number of non-admissible integers is a*b.")
    print("The only admissible k is k=a*b.")
    print("The final equation for the count of non-admissible integers is:")
    print(f"{a} * {b} = {result}")

# Case 2: Both dimensions are greater than 1.
else:
    # For this case, a known mathematical result provides the answer.
    # A theorem by Brezavšček, Oblak, and Šivic (2015) on orthogonal bases of matrices
    # states that for a, b > 1, an integer k is admissible if and only if k is not 1 and k is not (a*b - 1).
    # Therefore, there are exactly two non-admissible integers.
    result = 2
    non_admissible_2 = a * b - 1
    print(f"For a={a}, b={b}, since both dimensions are greater than 1, the number of non-admissible integers is 2.")
    print(f"The two non-admissible integers are 1 and (a*b - 1).")
    print(f"The first non-admissible integer is: 1")
    print(f"The equation for the second non-admissible integer is:")
    print(f"{a} * {b} - 1 = {non_admissible_2}")
