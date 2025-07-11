def display_final_equation():
  """
  This function prints the derived final expression for the synaptic efficacy dynamics.
  Note: The instruction "output each number in the final equation" is interpreted
  as displaying all components of the derived formula, as there are no numerical
  values given in the problem statement.
  """
  
  # Define unicode symbols for a cleaner representation of the equation
  tau_w = "τ_w"
  w_i_dot = "ẇᵢ"
  rho = "ρ"
  u_i = "uᵢ"
  v_i = "vᵢ"
  
  # Print the final equation using the symbols
  # The equation is: τ_w * ẇᵢ = ρ * uᵢ * (vᵢ / (1 + vᵢ))
  print("The final derived expression for the dynamics of synaptic efficacy is:")
  print(f"{tau_w} {w_i_dot} = {rho} {u_i} * ( {v_i} / (1 + {v_i}) )")

display_final_equation()