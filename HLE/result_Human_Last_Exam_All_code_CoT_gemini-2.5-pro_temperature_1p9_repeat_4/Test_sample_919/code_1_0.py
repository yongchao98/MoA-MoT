def solve_force_equation():
  """
  This function prints the derived symbolic equation for the force per unit area.
  """
  
  # Define the symbols used in the equation
  # f/area: Force per unit area vector
  # mu_0: Permeability of free space
  # mu: Permeability of the magnetic material
  # K_0: Amplitude of the surface current density
  # a: Spatial constant of the current
  # d: Thickness of the air gap
  # y: Spatial coordinate
  # i_x: Unit vector in the x-direction
  
  equation_string = "f_vec/area = (mu_0 / 2) * (K_0**2 * sin(a*y)**2) / (cosh(a*d) + (mu_0 / mu) * sinh(a*d))**2 * i_x_hat"
  
  print("The force per unit y-z area on the x = d interface is:")
  print(equation_string)

solve_force_equation()