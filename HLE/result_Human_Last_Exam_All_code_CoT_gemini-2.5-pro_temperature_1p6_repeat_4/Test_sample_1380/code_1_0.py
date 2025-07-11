def print_bachelier_discovery():
  """
  This function provides the answer to the user's historical question
  about Louis Bachelier's work.
  """
  discovery_name = "Brownian motion"
  description = (
      f"The specific physical process is {discovery_name}. "
      "Louis Bachelier, in his 1900 dissertation 'Théorie de la spéculation,' "
      "modeled the fluctuations of stock prices using a mathematical framework that was identical to the one "
      "describing the random, erratic movement of particles suspended in a fluid (e.g., pollen in water). "
      "This physical phenomenon is known as Brownian motion. Bachelier's work was the first to describe this "
      "process, now often called a Wiener process, with mathematical rigor, predating Albert Einstein's famous "
      "1905 paper on the subject."
  )
  print(description)

print_bachelier_discovery()