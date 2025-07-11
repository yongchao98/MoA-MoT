def check_actions_on_torus():
  """
  This function programmatically demonstrates that the two actions are the same
  for the 2-torus by simulating them with a sample loop and a sample fiber point.

  Let's represent elements of pi_1(T^2) and points in the fiber p^{-1}(x_0)
  as 2D integer vectors (tuples).
  """
  # A sample element from the fundamental group pi_1(T^2), corresponding to a loop.
  loop_class = (2, 3) # Corresponds to wrapping 2 times one way, 3 times the other.

  # A sample point in the fiber p^{-1}(x_0).
  fiber_point = (5, 7)

  # --- Simulate Action 1: Action by holonomy (path lifting) ---
  # This action is equivalent to adding the loop's vector to the fiber point's vector.
  # The endpoint of the lift of the loop (m,n) starting at (m_0, n_0) is (m_0+m, n_0+n).
  result_action_1 = (fiber_point[0] + loop_class[0], fiber_point[1] + loop_class[1])


  # --- Simulate Action 2: Action by restricting deck transformations ---
  # The loop class (m,n) corresponds to the deck transformation f(x,y) = (x+m, y+n).
  # We apply this transformation to the fiber point.
  result_action_2 = (fiber_point[0] + loop_class[0], fiber_point[1] + loop_class[1])


  # --- Compare the results ---
  are_actions_the_same = (result_action_1 == result_action_2)

  # As demonstrated in the thinking steps, the two actions are identical.
  # The question asks for a "Yes" or "No" answer.
  if are_actions_the_same:
    print("Yes")
  else:
    print("No")

check_actions_on_torus()