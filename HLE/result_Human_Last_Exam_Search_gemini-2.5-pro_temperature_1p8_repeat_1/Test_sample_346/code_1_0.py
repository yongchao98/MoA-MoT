from dovekie import S, K, I
import math

# Inside-out construction to avoid parenthesis errors
ssksk = S(S(K(S))(K))
sii = S(I)(I)
f = ssksk(I)
kf = K(f)
fg = ssksk(sii)
si_si_kf = S(I)(S(I)(kf))
si_si_kf_fg = si_si_kf(fg)

# The full expression. This is based on re-analyzing the structure.
# S(I) applies its second argument to the result of its first argument applied to the second.
# S(I)(A)(B) -> B(A(B))
# Here A = S(I)(S(I)(K f)) and B = f g
# So expr = fg ( (S(I)(S(I)(K f))) (fg) )
# The online interpreter shows S(I)(A)(B) reduces to A B (I B) -> A B B
# No, S x y z = x z (y z), so S(I) a b = I b (a b) = b (a b). This is correct.

# I suspect the original prompt's expression is written in a slightly non-standard way.
# Let's try to interpret the application order as left-associative.
expr_text = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"

# A will be S(I)
# B will be S(I)(S(I)(K(S(S(K(S))(K))(I)))))
# C will be S(S(K(S))(K))(S(I)(I))
# So it becomes A(B)(C). This doesn't match S's arity.

# Let's try again with the `dovekie` code, paying extreme attention to the parenthesis.