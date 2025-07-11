import math

# The problem asks for the exact number of F_q-rational maximal tori
# of a reductive group G of type E_8 over the finite field F_q.

# The number of F_q-rational maximal tori of a reductive group G over F_q is equal to
# the number of F-conjugacy classes in its Weyl group W, where F is the Frobenius endomorphism.
# For a group of type E_8, any such group over a finite field is split. This implies
# that the Frobenius map F acts trivially on the Weyl group W.
# Therefore, the F-conjugacy classes are simply the ordinary conjugacy classes of W.
# The problem is thus equivalent to finding the number of conjugacy classes of the Weyl group W(E_8).

# The number of conjugacy classes of W(E_8) is a known result from the theory of Weyl groups.
# A direct computation is computationally very expensive due to the large order of W(E_8),
# which is |W(E_8)| = 2^14 * 3^5 * 5^2 * 7 = 696,729,600.
# The number of its conjugacy classes has been established to be 112.

num_tori = 112

print("Let G be a reductive group of type E_8 over the finite field F_q.")
print("The number of F_q-rational maximal tori of G corresponds to the number of conjugacy classes in its Weyl group, W(E_8).")
print("\nFinal equation:")
# The following line outputs the numbers in the final equation format as requested.
print(f"Number of F_q-rational maximal tori = Number of conjugacy classes in W(E_8) = {num_tori}")
