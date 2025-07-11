def describe_process():
  """
  This function describes the biological process shown in the provided image.
  """
  process_name = "Reorientation of tip growth"
  description = (
      "The image sequence (L, M, N) illustrates the process of reorientation of tip growth in a filamentous organism.\n"
      "1. In panel L, the filament is growing straight.\n"
      "2. In panel M, a new growth point initiates as a bulge on the side of a cell (indicated by the arrowhead).\n"
      "3. In panel N, this new point has developed, causing the filament to bend and change its overall direction of growth."
  )
  
  print(f"Process: {process_name}\n")
  print("Description:")
  print(description)

describe_process()