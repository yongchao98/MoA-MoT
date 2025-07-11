import math

def gcd(a, b):
  """Computes the greatest common divisor of a and b."""
  return math.gcd(a, b)

def main():
  """
  Calculates the dimension of the cohomology group H^2(G,M).
  """
  dim_M = 128
  order_a = float('inf')
  order_b = float('inf')
  
  # The action of a and b is a cyclic permutation of order 128.
  cycle_len = 128
  
  # Step 1: Calculate dimensions of cohomology groups for A, B, and H.
  # For a cyclic group <t> acting as sigma^p on M=Q^n where sigma is an n-cycle,
  # the dimension of the invariant subspace M^<t> is gcd(p, n).
  # For Z, dim H^k(Z, M) = dim(M^<t>) for k=0,1. And H^2 = H^0.
  
  # For A = <a>, the action is sigma^1.
  p_A = 1
  dim_Hk_A = gcd(p_A, cycle_len)
  
  # For B = <b>, the action is sigma^1.
  p_B = 1
  dim_Hk_B = gcd(p_B, cycle_len)

  # For H = <a^8> = <b^8>, the action is sigma^8.
  p_H = 8
  dim_Hk_H = gcd(p_H, cycle_len)

  print(f"Plan step 1: Calculate dimensions of cohomology groups.")
  print(f"dim H^k(A, M) for k=0,1,2 is gcd({p_A}, {cycle_len}) = {dim_Hk_A}")
  print(f"dim H^k(B, M) for k=0,1,2 is gcd({p_B}, {cycle_len}) = {dim_Hk_B}")
  print(f"dim H^k(H, M) for k=0,1,2 is gcd({p_H}, {cycle_len}) = {dim_Hk_H}")
  print("-" * 20)
  
  # Step 2: Analyze f1 to find dim(coker(f1)).
  # f1: H^1(A,M) + H^1(B,M) -> H^1(H,M)
  # dim(Domain(f1)) = dim_Hk_A + dim_Hk_B = 1 + 1 = 2
  # dim(Codomain(f1)) = dim_Hk_H = 8
  # The map f1 sends (u,v) to res_A(u) - res_B(v).
  # Since the actions of a and b are identical, the restriction maps are identical.
  # So Image(f1) = Image(res_A).
  # The restriction map R: H^1(A,M) -> H^1(H,M) is non-zero.
  # Since dim H^1(A,M) = 1, the image of R must have dimension 1.
  dim_Im_f1 = 1
  dim_coker_f1 = dim_Hk_H - dim_Im_f1

  print(f"Plan step 2: Analyze the map f1.")
  print(f"The map f1 has domain dimension {dim_Hk_A + dim_Hk_B} and codomain dimension {dim_Hk_H}.")
  print(f"The image of f1 has dimension {dim_Im_f1}.")
  print(f"dim(coker(f1)) = dim H^1(H,M) - dim(Im(f1)) = {dim_Hk_H} - {dim_Im_f1} = {dim_coker_f1}")
  print("-" * 20)

  # Step 3: Analyze f2 to find dim(ker(f2)).
  # f2: H^2(A,M) + H^2(B,M) -> H^2(H,M)
  # H^2(K,M) = H^0(K,M) = M^K for K=Z.
  # f2 maps (u,v) to res(u) - res(v).
  # The restriction map res: M^A -> M^H is just an inclusion map.
  # So f2(u,v) = u - v.
  # The kernel of f2 is the set of pairs (u,u) where u is in H^2(A,M).
  # So dim(ker(f2)) = dim(H^2(A,M)).
  dim_ker_f2 = dim_Hk_A

  print(f"Plan step 3: Analyze the map f2.")
  print(f"The map f2 has domain dimension {dim_Hk_A + dim_Hk_B} and codomain dimension {dim_Hk_H}.")
  print(f"dim(ker(f2)) is equal to dim H^2(A,M) which is {dim_ker_f2}.")
  print("-" * 20)

  # Step 4: Combine the results.
  # dim H^2(G,M) = dim(coker(f1)) + dim(ker(f2))
  dim_H2_G_M = dim_coker_f1 + dim_ker_f2
  
  print(f"Plan step 4: Final calculation.")
  print(f"From the Mayer-Vietoris sequence, we have the equation:")
  print(f"dim H^2(G,M) = dim(coker(f1)) + dim(ker(f2))")
  print(f"dim H^2(G,M) = ({dim_Hk_H} - {dim_Im_f1}) + {dim_ker_f2}")
  print(f"dim H^2(G,M) = {dim_coker_f1} + {dim_ker_f2} = {dim_H2_G_M}")

if __name__ == '__main__':
  main()
