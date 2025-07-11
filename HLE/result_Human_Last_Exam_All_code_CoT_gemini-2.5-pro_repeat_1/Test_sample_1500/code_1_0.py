# A demonstration for the counterexample in part (b).
# We consider the case G = SU(n) with n=3.
n = 3

# The question asks if the second Betti number b_2(O_lambda) is always n-1.
# Let's calculate the value predicted by this statement.
predicted_b2 = n - 1

# Now, consider a specific coadjoint orbit O_lambda for SU(3).
# If lambda is a singular element (on the boundary of the Weyl alcove),
# for example, the highest weight of the standard representation,
# the corresponding orbit is the complex projective plane CP^2.
# The second Betti number of CP^2 is a well-known topological invariant.
actual_b2_for_CP2 = 1

# Now we check if the prediction holds for this case.
is_statement_true_for_this_case = (actual_b2_for_CP2 == predicted_b2)

print("--- Analysis for Question (b) ---")
print(f"For G = SU(n) with n = {n}:")
print(f"The statement in (b) suggests b_2(O_lambda) should always be n - 1 = {predicted_b2}.")
print("However, we can choose a lambda on the boundary of the Weyl alcove.")
print(f"For n = {n}, such an orbit can be the complex projective space CP^2.")
print(f"The actual second Betti number for CP^2 is b_2(CP^2) = {actual_b2_for_CP2}.")
print(f"Comparing the values: Is {actual_b2_for_CP2} == {predicted_b2}? {is_statement_true_for_this_case}.")
print("Since the values do not match, the statement in (b) is false.")
print("\nThe final equation is:")
print(f"b_2(CP^2) = {actual_b2_for_CP2}")
print(f"n - 1 = {predicted_b2}")
print(f"The comparison shows {actual_b2_for_CP2} != {predicted_b2}.")